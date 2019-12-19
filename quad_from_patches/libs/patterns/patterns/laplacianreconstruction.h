#ifndef LAPLACIANRECONSTRUCTION_H
#define LAPLACIANRECONSTRUCTION_H

#include <Eigen/Sparse>
#include <memory>
#include <unordered_map>
#include <vcg/complex/algorithms/mesh_to_matrix.h>
#include <vcg/complex/algorithms/implicit_smooth.h>
#include <vcg/complex/algorithms/update/quality.h>

template<class MeshType>
class LaplacianReconstruction
{
public:
    typedef typename MeshType::VertexType     VertexType;
    typedef typename MeshType::VertexPointer  VertexPointer;
    typedef typename MeshType::VertexIterator VertexIterator;
    typedef typename MeshType::FaceType       FaceType;
    typedef typename MeshType::FacePointer    FacePointer;
    typedef typename MeshType::FaceIterator   FaceIterator;
    typedef typename MeshType::ScalarType     ScalarType;
    typedef typename MeshType::CoordType      CoordType;
    typedef Eigen::Matrix<ScalarType,  3, 1> Value ;
    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 3> Vector;
    typedef Eigen::SparseMatrix <ScalarType> Matrix;
    int numberConstraints;

    // members
    ScalarType laplace_constraintWeight;
    Matrix L;
    Vector b;
    std::shared_ptr<Eigen::SimplicialCholesky<Matrix>> solver;
    MeshType * ptrmesh;
    std::vector<int> indexes_constraints;

    LaplacianReconstruction(){
        laplace_constraintWeight=10000;
        numberConstraints=0;
        indexes_constraints.clear();
        ptrmesh=NULL;
    }
    LaplacianReconstruction(MeshType & mesh){
        vcg::tri::RequireCompactness(mesh);
        ptrmesh=&mesh;
        laplace_constraintWeight=10000;
        indexes_constraints.clear();
        if(setConstraintsFromSelectedVertices())
            constructL();
    }
    void setTopology(MeshType &mesh){
        vcg::tri::RequireCompactness(mesh);
        ptrmesh=&mesh;

        if(setConstraintsFromSelectedVertices())
            constructL();
    }
    bool setConstraintsFromSelectedVertices(){
        //int numberSelected=0;
        indexes_constraints.clear();
        int nv=ptrmesh->VN();
        for(VertexIterator vi=ptrmesh->vert.begin();vi!=ptrmesh->vert.end();vi++){
            if(vi->IsS()){
                //numberSelected++;
                indexes_constraints.push_back(vcg::tri::Index(*ptrmesh,&*vi));
            }
        }
        if(indexes_constraints.empty())
            return false;
        else{
            numberConstraints=indexes_constraints.size();
            // set right hand side of the system
            //b.resize(L.rows(),N);
            b.setZero(nv+numberConstraints,3);
            for(int i=0;i<indexes_constraints.size();i++){
                    CoordType point=ptrmesh->vert[indexes_constraints[i]].P();
                    b.row(nv+i)<<laplace_constraintWeight*point.X() , laplace_constraintWeight*point.Y() , laplace_constraintWeight*point.Z();
            }
            return true;
        }
    }
    void constructL(){
        //get the entries for laplacian matrix
        std::vector<std::pair<int,int> > IndexValence;
        std::unordered_map<int,set<int>> startVV;
        //vcg::tri::MeshToMatrix<MeshType>::GetLaplacianMatrix(*ptrmesh,IndexL,ValuesL,false,1,false);
        UpdateTopology<MeshType>::FaceFace(*ptrmesh);
        //UpdateQuality<MeshType>::VertexValence(*ptrmesh);
        UpdateFlags<MeshType>::VertexClearV(*ptrmesh);
        UpdateFlags<MeshType>::FaceBorderFromFF(*ptrmesh);
        for(FaceIterator fi=ptrmesh->face.begin();fi!=ptrmesh->face.end();fi++){
            if(!fi->IsD()){
                for(int i=0;i<fi->VN();i++){
                    if(!fi->V(i)->IsV()){
                        vcg::face::Pos<FaceType> pos(&*fi,i,fi->V(i));
                        vector<vcg::face::Pos<FaceType>> posVec;
                        set<int> areedges;
                        int index=vcg::tri::Index(*ptrmesh,fi->V(i));
                        int valence=0;
                        bool foundBorder=false;

                          vcg::face::Pos<FaceType> curPos=pos;
                          do
                          {
                            assert(curPos.IsManifold());
                            if(curPos.IsBorder() && !foundBorder) {
                              foundBorder=true;
                            }
                            posVec.push_back(curPos);
                            areedges.insert(vcg::tri::Index(*ptrmesh,curPos.VFlip()));
                            curPos.FlipF();
                            curPos.FlipE();
                          } while(curPos!=pos);
                        valence=posVec.size();
                        if(foundBorder){
                            valence=posVec.size()/2+1;
                        }
                        assert(valence==areedges.size());
                        IndexValence.push_back(make_pair(index,valence));
                        startVV.emplace(index,areedges);
                        fi->V(i)->SetV();
                    }
                }
            }
        }


        int nv=IndexValence.size();
        int m=numberConstraints+nv;
        std::vector<Eigen::Triplet<ScalarType> > IJV;
        IJV.reserve(m);
        for(size_t i= 0;i<IndexValence.size();i++)
        {
            int index=IndexValence[i].first;
            int di=IndexValence[i].second;
            set<int> edges=startVV[index];
            assert(index<m);
            for(set<int>::iterator sit=edges.begin();sit!=edges.end();sit++){
                assert(*sit<m);
                IJV.push_back(Eigen::Triplet<ScalarType>(index,*sit,-1));
            }
            IJV.push_back(Eigen::Triplet<ScalarType>(index,index,di));
        }
        for(int i=0;i<numberConstraints;i++){
            IJV.push_back(Eigen::Triplet<ScalarType>(nv+i,indexes_constraints[i],laplace_constraintWeight));
        }
        L.resize(m,nv);
        L.setFromTriplets(IJV.begin(),IJV.end());
        //cout<<"l size "<< L.rows()<<" "<<L.cols()<< endl;
        //cout<<L<<endl;
        //cout<<"b"<<endl;
        //cout<<b<<endl;
        solver.reset(new Eigen::SimplicialCholesky<Matrix>(L.transpose() * L));
    }

    void laplace_solve() {
        Vector x = solver->solve(L.transpose()*b);
        int nv=ptrmesh->VN();
        // copy result to Data::value
        for (int i = 0; i < nv; ++i) {
            if (!ptrmesh->vert[i].IsS()){
                ptrmesh->vert[i].P()=CoordType(x(i,0),x(i,1),x(i,2));
            }
        }
    }
};

#endif // LAPLACIANRECONSTRUCTION_H
