#ifndef PATCHG_H
#define PATCHG_H
#include "laplacianreconstruction.h"
#include "ktmethod/patchgen/decl.h"
#include "ktmethod/patchgen/generate_subtopology.h"
#include <wrap/io_trimesh/export_obj.h>
#include "myutils.h"
//#include "ktmethod/patchgen/generate_topology.h"

template <class MeshType,class PatchParam>
class PatchG
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
    MeshType mesh;
    LaplacianReconstruction<MeshType> reconstructor;
    bool finish;
    int numberSides;
    vector<int> paddingValues;
    vector<int> borders;
    vector<int> corners;
    //for debug
    vector<CoordType> points;
    PatchG(){
        mesh.Clear();
        finish=false;
        numberSides=0;
        paddingValues.clear();
    }
    void generate_topology(const PatchParam& param) {
        int num_sides = param.get_num_sides();
        paddingValues.clear();
        finish=false;
        mesh.Clear();
        numberSides=num_sides;
//        std::cout << "nsides =" << num_sides <<endl;
//        std::cout << "pattern=" << param.pattern_id;
        for(int i=0;i<num_sides;i++)
            paddingValues.push_back(param.p[i]);

        std::string param_str = patchgen::get_param_str(num_sides, param.pattern_id, param);
//        if (!param_str.empty())
//            std::cout << ", param=" << param_str;
//        std::cout << "\n";

        // generate topology
        patchgen::generate_subtopology(num_sides, param.pattern_id, param, mesh);
        /*determine_geometry(param.l);
        int Mask=0;
        string aux="meshbeforepadding"+std::to_string(mesh.vert.size()+mesh.face.size())+".obj";
        vcg::tri::io::ExporterOBJ<MeshType>::Save(mesh,aux.c_str(), Mask);
        */
        //UpdateNormal<MeshType>::PerPolygonalFaceNormalized(mesh);
        addPading();

        // reordering corner index, flip faces
        /*for (auto v : patch.vertices()) {
            int& corner_index = patch.data(v).patchgen.corner_index;
            if (corner_index == -1) continue;
            corner_index = param.permutation[corner_index];
            if (param.permutation.is_flipped()) corner_index = (corner_index + 1) % num_sides;      // Be careful!
        }
        if (param.permutation.is_flipped())
            kt84::flip_faces(patch);*/
    }

    void generate_topology(const Eigen::VectorXi& l, PatchParam& param) {
        int num_sides = l.size();
//        std::cout << "input=[";
//        for (int i = 0; i < num_sides; ++i)
//            std::cout << " " << l[i];
//        std::cout << " ]\n";

        param = patchgen::get_default_parameter(l);
//        std::cout << "\n";

        generate_topology(param);
    }

    void determine_geometry(const Eigen::VectorXi& l) {
        // fix boundary vertices
        UpdateTopology<MeshType>::FaceFace(mesh);
        UpdateFlags<MeshType>::FaceBorderFromFF(mesh);

        int num_sides = numberSides;

        vcg::face::Pos<FaceType> start;
        start.SetNull();
        for(typename MeshType::FaceIterator fi=mesh.face.begin();fi!=mesh.face.end();fi++){
            if(!fi->IsD()){
                for(int k=0;k<fi->VN();k++){
                    if(vcg::face::IsBorder(*fi,k)){
                         start.Set(&*fi,k,fi->V(k));
                         break;
                    }
                }
            }
        }
        vector<VertexPointer> frontiers;
        frontiers.push_back(start.V());
        do{
            do{
                start.FlipE();
                start.FlipF();
            }while(!start.IsBorder());
            start.FlipV();
            frontiers.push_back(start.V());
        }while(start.V()!=frontiers[0]);
        frontiers.pop_back();

        Eigen::VectorXi current(numberSides);
        for(int i=0;i<numberSides;i++)
            current[i]=get_SideK(mesh,i).size()-1;
//        cout<<"current version "<<current.transpose()<<endl;
        points.clear();
        if(num_sides>2){
            double side_angle = 2 * 3.1416 / num_sides;
            int t=0;
            for(int i=0;i<frontiers.size();i++){
                if(frontiers[i]->IsS()){
                    double theta=side_angle*t;
                    frontiers[i]->P()=CoordType(cos(theta),sin(theta),0.0);
//                    cout<<"corner index "<<vcg::tri::Index(mesh,frontiers[i])<<endl;
                    this->corners.push_back(vcg::tri::Index(mesh,frontiers[i]));
                    t++;
                }
            }

            for(int i=0;i<num_sides;i++){
                vector<int> sideIndexes=get_SideK(mesh,i);
                int length=sideIndexes.size()-1;
                assert(mesh.vert[sideIndexes[0]].IsS());
                assert(mesh.vert[sideIndexes.back()].IsS());
                assert(length>0);
                //cout<<"side indexes "<<sideIndexes[0]<<" "<<sideIndexes.back()<<endl;
                CoordType e0=mesh.vert[sideIndexes[0]].P();
                CoordType e1=mesh.vert[sideIndexes.back()].P();

                //cout<<"corner point "<<e0.X()<<" "<<e0.Y()<<" "<<e0.Z()<<endl;
                //cout<<"corner point "<<e1.X()<<" "<<e1.Y()<<" "<<e1.Z()<<endl;
                for(int j=1;j<sideIndexes.size()-1;j++){
                    float m=float(j)/length;
                    float n=float(length-j)/length;
                   // cout<<"m "<<m<<" n "<<n<<endl;
                    //e0.Scale(n,n,n);
                    //e1.Scale(m,m,m);
                    mesh.vert[sideIndexes[j]].P()=e0*n+e1*m;
                    //points.push_back(mesh.vert[sideIndexes[j]].P());
                    //cout<<"point "<<mesh.vert[sideIndexes[j]].P().X()<<" "<<mesh.vert[sideIndexes[j]].P().Y()<<" "<<mesh.vert[sideIndexes[j]].P().Z()<<endl;
                    //cout<<"point "<<mesh.vert[sideIndexes[j]].P().X()<<" "<<mesh.vert[sideIndexes[j]].P().Y()<<" "<<mesh.vert[sideIndexes[j]].P().Z()<<endl;
                    mesh.vert[sideIndexes[j]].SetS();
                }
            }

        }
        else if(num_sides>1){

            num_sides =frontiers.size();
            double side_angle = 2 * 3.1416 / num_sides;
            for(int i=0;i<num_sides;i++){
                double theta = side_angle * i;
                CoordType vdata(cos(theta),sin(theta),0);
                frontiers[i]->P()=vdata;
                frontiers[i]->SetS();
            }
        }

        // solve

        reconstructor.setTopology(mesh);
        reconstructor.laplace_solve();

        finish=true;
        UpdateNormal<MeshType>::PerPolygonalFaceNormalized(mesh);
        // correct orientation some faces
        for(FaceIterator fi=mesh.face.begin();fi!=mesh.face.end();fi++){
            if(!fi->IsD()){
                if(fi->N().Z()<0){
                    std::vector<VertexType*> faceVert;
                    for (int i=0;i<fi->VN();i++)
                        faceVert.push_back(fi->V(i));
                    std::reverse(faceVert.begin(),faceVert.end());
                    for (int i=0;i<fi->VN();i++)
                        fi->V(i)=faceVert[i];
                }
            }
        }

        UpdateTopology<MeshType>::FaceFace(mesh);
        UpdateFlags<MeshType>::FaceBorderFromFF(mesh);
        UpdateNormal<MeshType>::PerPolygonalFaceNormalized(mesh);

        //restore corners
        UpdateFlags<MeshType>::VertexClearS(mesh);
        for(int i=0;i<corners.size();i++)
            mesh.vert[corners[i]].SetS();

        start.SetNull();
        for(typename MeshType::FaceIterator fi=mesh.face.begin();fi!=mesh.face.end();fi++){
            if(!fi->IsD()){
                assert(fi->N().Z() > 0);
                for(int k=0;k<fi->VN();k++){
                    if(vcg::face::IsBorder(*fi,k)){
                         start.Set(&*fi,k,fi->V(k));
                         break;
                    }
                }
            }
        }

        vector<VertexPointer> bordersPointers;
        bordersPointers.push_back(start.V());
        do{
            do{
                start.FlipE();
                start.FlipF();
            }while(!start.IsBorder());
            start.FlipV();
            bordersPointers.push_back(start.V());
        }while(start.V()!=bordersPointers[0]);
        bordersPointers.pop_back();

        this->corners.clear();
        this->borders.clear();
        for (VertexPointer vp : bordersPointers) {
            if (vp->IsS()) {
                this->corners.push_back(vcg::tri::Index(mesh,vp));
            }
            this->borders.push_back(vcg::tri::Index(mesh,vp));
        }

        //Clean<MeshType>::OrientCoherentlyMesh(mesh,isoriented,isorientable);
        //Clean<MeshType>::FlipMesh(mesh);
    }

    void determine_geometry(const vector<vector<typename MeshType::CoordType>> &l) {
        assert(numberSides==l.size());
        // fix boundary vertices
        UpdateTopology<MeshType>::FaceFace(mesh);
        UpdateFlags<MeshType>::FaceBorderFromFF(mesh);


        auto side_indexL=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("LeftSide"));
        auto side_indexR=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("RightSide"));

       // reordering the input. This because the pattern may have been permuted
        vector<vector<int>> sideIndexes(numberSides);
        vector<int> sideIndexes_l(numberSides);
        vector<int> l_l(numberSides);
        for(int i=0;i<numberSides;i++){
            sideIndexes[i]=getSideK(i);
            sideIndexes_l[i]=sideIndexes[i].size()-1;
            l_l[i]=l[i].size()-1;
        }
        /*determine_geometry(Eigen::VectorXi::Zero(5));
        int Mask=0;
        string aux="../exported/patchmesh2D"+std::to_string(numberSides)+"-"+std::to_string(mesh.vert.size())+".obj";
        vcg::tri::io::ExporterOBJ<MeshType>::Save(mesh,aux.c_str(), Mask);*/

        //UpdateFlags<MeshType>::VertexClearS(mesh);
        std::vector<std::pair<int,int>> correspondence;
        if(!utility::ispermuted(l_l,sideIndexes_l,correspondence)){
            finish=false;
//            cout<<"INCONSISTENTE PATCHES 2D"<<endl;
            return;
        }
        /*vector<int> corners=getCorners();
        vector<typename MeshType::CoordType> positions;
        positions.push_back(l[0][0]);
        positions.push_back(l[0].back());
        for(int i=1;i<numberSides;i++){
            if(positions.back()==l[i][0])
                positions.push_back(l[i].back());
            else
                positions.push_back(l[i][0]);
        }
        positions.pop_back();
        assert(positions.size()==numberSides);*/
        if(sideIndexes[correspondence[0].second].back()!=sideIndexes[correspondence[1].second][0])
            for(int i=0;i<sideIndexes.size();i++){
                vector<int> newsind=sideIndexes[i];
                std::reverse(newsind.begin(),newsind.end());
                sideIndexes[i]=newsind;
            }

        for(int i=0;i<numberSides;i++){
            vector<int> positions=sideIndexes[correspondence[i].second];            
            for(int j=0;j<positions.size();j++){
                //cout<<"mesh index "<<positions[j]<<", index "<<j<<" associated [ "<<l[i][j].X()<<" "<<l[i][j].Y()<<endl;
                mesh.vert[positions[j]].P()=l[i][j];
                mesh.vert[positions[j]].SetS();
            }
        }

        // solve

        reconstructor.setTopology(mesh);
        reconstructor.laplace_solve();

        // Correcting points that remain outside the boundary of the patch (due to the large difference between the number of divisions on one side of the patch and the others)
        // we attempt to walk by the boundary and check tha adjacent vertices

        // converting the boundary to a vector of points
        /*vector<CoordType> boundary;
        for(size_t i=0;i<l.size();i++){
            vector<CoordType> temp=l[i];
            temp.pop_back();
            boundary.insert(boundary.end(),temp.begin(),temp.end());
        }
        face::Pos<FaceType> startPosBoundary;
        startPosBoundary.SetNull();
        for(auto fi=mesh.face.begin();fi!=mesh.face.end();++fi){
            if(!(*fi).IsD()){
                for(int i=0;i<fi->VN();++i)
                {
                    if(vcg::face::IsBorder(*fi,i)){
                        startPosBoundary.Set(&*fi,i,fi->V(i)); // start Pos edge
                        break;
                    }
                }
                if(!startPosBoundary.IsNull())
                    break;
            }
        }
        face::Pos<FaceType> copy=startPosBoundary;
        do{
            face::Pos<FaceType> next=copy;
            do{
                next.FlipE();
                next.FlipF();
            }while(!next.IsBorder());
            copy.FlipE();
            if(!utility::pointIsInsidePoly(boundary,copy.VFlip()->P())){
               cout<<"a point detected outside the polygon"<<endl;
               //computing reflexion over copy and next
               copy.FlipE();
               CoordType P=copy.VFlip()->P();
               CoordType reflexion1,reflexion2; //  actually, we must be check the valence equal 2, but for while ...
               utility::computeReflexion(reflexion1,P,copy.V()->P(),copy.VFlip()->P());
               utility::computeReflexion(reflexion1,P,copy.V()->P(),copy.VFlip()->P());
            }
            next.FlipV();
            copy=next;
        }while(copy!=startPosBoundary);

        finish=true;

        UpdateNormal<MeshType>::PerPolygonalFaceNormalized(mesh);
        // correct orientation some faces
        for(FaceIterator fi=mesh.face.begin();fi!=mesh.face.end();fi++){
            if(!fi->IsD()){
                if(fi->N().Z()<0){
                    VertexPointer temp=fi->V(3);
                    fi->V(3)=fi->V(1);
                    fi->V(1)=temp;
                }
            }
        }
        UpdateTopology<MeshType>::FaceFace(mesh);
        UpdateNormal<MeshType>::PerPolygonalFaceNormalized(mesh);*/
    }
    void create_and_process(const vector<vector<typename MeshType::CoordType>> &borderpatch){
         Eigen::VectorXi l(borderpatch.size());
         //cout<<"creating patch with "<<borderpatch.size()<<" sides "<<endl;
         for(int i=0;i<borderpatch.size();i++){
             l[i]=borderpatch[i].size()-1;
             //cout<<"side "<<l[i]<<endl;
         }
//         if(l.sum()%2!=0)
//             cout<<"falhou"<<endl;
         assert(l.sum()%2==0);
         vector<vector<typename MeshType::CoordType>> positions=borderpatch;
         /*if(l.sum()%2!=0){
             l[0]+=1;
             vector<typename MeshType::CoordType> temp=positions[0];
             temp.insert(temp.begin(),positions[0][0]);
             positions[0]=temp;
         }*/

         PatchParam param = patchgen::get_default_parameter(l);
         generate_topology(param);
         determine_geometry(positions);
    }
    vector<int> getSideK(int side){
        assert(side<numberSides);
        vcg::tri::RequireCompactness(mesh);
        vector<int> result;
        auto side_indexL=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("LeftSide"));
        auto side_indexR=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("RightSide"));
        bool founded=false;
        vcg::face::Pos<FaceType> start;
//        cout<<"side "<<side<<endl;
        for(FaceIterator fi=mesh.face.begin();fi!=mesh.face.end();fi++){
            for(int k=0;k<fi->VN();k++){
                if(fi->V(k)->IsS()){
//                    cout<<"corner "<<tri::Index(mesh,fi->V(k))<<endl;
//                    cout<<"sideL "<<side_indexL[fi->V(k)]<<endl;
//                    cout<<"sideR "<<side_indexR[fi->V(k)]<<endl;
                    if(side_indexR[fi->V(k)]==side){
//                        cout<<"corner "<<tri::Index(mesh,fi->V(k))<<" is right"<<endl;
//                        cout<<"searching extreme"<<endl;
                        start.Set(&*fi,k,fi->V(k));
                        int count=0;
                        while(side_indexL[start.VFlip()]!=side && count<10){
                            start.NextE();
//                            cout<<"test "<<tri::Index(mesh,start.VFlip())<<endl;
//                            cout<<" sideL "<<side_indexL[start.VFlip()]<<endl;
//                            cout<<" sideR "<<side_indexR[start.VFlip()]<<endl;
                            count++;
                        }
                        start.FlipV();
//                        cout<<"extreme "<<tri::Index(mesh,fi->V(k))<<endl;
                        founded=true;
                        break;
                    }
                }
            }
            if(founded)
                break;
        }
        result.push_back(vcg::tri::Index(mesh,start.VFlip()));
        vcg::face::Pos<FaceType> copy=start;
        do{
            result.push_back(vcg::tri::Index(mesh,start.V()));            
            do{
                start.NextE();
            }while(side_indexL[start.VFlip()]!=side && side_indexR[start.VFlip()]!=side);
            if(start==copy)
                break;
            start.FlipV();
            if(side_indexL[start.V()]!=side_indexR[start.V()]){
                result.push_back(vcg::tri::Index(mesh,start.V()));
                break;
            }
        }while(side_indexL[start.V()]==side_indexR[start.V()]);

//        cout<<"side "<<side<<endl;
//        for(int i=0;i<result.size();i++)
//            cout<<" "<<result[i];
//        cout<<endl;
        return result;
    }
    void addPading(){
        assert(!paddingValues.empty());
        for(int i=0;i<paddingValues.size();i++){
//            cout<<"p["<<i<<"]= "<<paddingValues[i]<<endl;
            for (int j = 0; j < paddingValues[i];j++)
                addQuadStrip_on_sideK(mesh,(numberSides-i)%numberSides);
        }
    }
};

#endif // PATCHG_H
