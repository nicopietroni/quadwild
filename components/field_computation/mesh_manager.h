/***************************************************************************/
/* Copyright(C) 2021


The authors of

Reliable Feature-Line Driven Quad-Remeshing
Siggraph 2021


 All rights reserved.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
****************************************************************************/

#ifndef MESH_MANAGER
#define MESH_MANAGER

//#define MINDOT -0.99

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/simplex/face/topology.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
//#include <wrap/gl/trimesh.h>
#include "AutoRemesher.h"
#include <wrap/io_trimesh/export_field.h>
#include <iostream>
#include <fstream>
#include <vcg/complex/algorithms/attribute_seam.h>
#include <vcg/complex/algorithms/crease_cut.h>
#include "mesh_field_smoother.h"
#include <vcg/complex/algorithms/polygonal_algorithms.h>

#include "fields/field_smoother.h"

// Basic subdivision class
template <class FaceType>
struct SplitLev : public   std::unary_function<vcg::face::Pos<FaceType> ,typename FaceType::CoordType >
{
    typedef typename FaceType::CoordType CoordType;
    typedef typename FaceType::VertexType VertexType;
    typedef typename FaceType::ScalarType ScalarType;

    typedef std::pair<CoordType,CoordType> CoordPair;
    std::map<CoordPair,CoordType> *SplitOps;

    void operator()(VertexType &nv,vcg::face::Pos<FaceType>  ep)
    {
        VertexType* v0=ep.f->V0(ep.z);
        VertexType* v1=ep.f->V1(ep.z);

        assert(v0!=v1);

        CoordType Pos0=v0->P();
        CoordType Pos1=v1->P();

        CoordPair CoordK(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
        assert(SplitOps->count(CoordK)>0);
        nv.P()=(*SplitOps)[CoordK];
    }

    vcg::TexCoord2<ScalarType> WedgeInterp(vcg::TexCoord2<ScalarType> &t0,
                                           vcg::TexCoord2<ScalarType> &t1)
    {
        return (vcg::TexCoord2<ScalarType>(0,0));
    }

    SplitLev(std::map<CoordPair,CoordType> *_SplitOps){SplitOps=_SplitOps;}
};

template <class FaceType>
class EdgePred
{
    typedef typename FaceType::CoordType CoordType;
    typedef typename FaceType::VertexType VertexType;
    typedef std::pair<CoordType,CoordType> CoordPair;
    std::map<CoordPair,CoordType> *SplitOps;

public:

    bool operator()(vcg::face::Pos<FaceType> ep) const
    {
        VertexType* v0=ep.f->V0(ep.z);
        VertexType* v1=ep.f->V1(ep.z);

        assert(v0!=v1);

        CoordType Pos0=v0->P();
        CoordType Pos1=v1->P();

        CoordPair CoordK(std::min(Pos0,Pos1),std::max(Pos0,Pos1));

        return (SplitOps->count(CoordK)>0);
    }

    EdgePred(std::map<CoordPair,CoordType> *_SplitOps){SplitOps=_SplitOps;}
};


template <class MeshType>
class MeshPrepocess
{
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType  CoordType;

    typedef typename MeshType::VertexType    VertexType;
    typedef typename MeshType::VertexPointer VertexPointer;

    typedef typename MeshType::FaceType    FaceType;
    typedef typename MeshType::FacePointer FacePointer;

    typedef std::pair<CoordType,CoordType> CoordPair;

    static bool IsFold(const FaceType &f,
                       const size_t &IndexE,
                       ScalarType MinDot=-0.99)
    {
        if (vcg::face::IsBorder(f,IndexE))return false;
        FaceType *fOpp=f.cFFp(IndexE);
        //int IOpp=f.cFFi(IndexE);
        CoordType N0=f.cN();
        CoordType N1=fOpp->N();
        N0.Normalize();
        N1.Normalize();
        if ((N0*N1)>MinDot)return false;
        return true;
    }

    typedef SplitLev<FaceType> SplitLevType;
    typedef EdgePred<FaceType> EdgePredType;

    static bool SplitFolds(MeshType &mesh,
                           ScalarType MinDot=-0.99,
                           bool debugmsg=false)
    {
        mesh.InitFeatureCoordsTable();

        //then save the edges to be splitted
        std::map<CoordPair,CoordType> ToBeSplitted;
        for (size_t i=0;i<mesh.face.size();i++)
        {
            //find the number of edges
            for (size_t j=0;j<3;j++)
            {
                if (!IsFold(mesh.face[i],j,MinDot))continue;
                int VIndex0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                int VIndex1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                CoordType P0=mesh.vert[VIndex0].P();
                CoordType P1=mesh.vert[VIndex1].P();
                std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
                ToBeSplitted[Key]=(P0+P1)/2;
            }
        }
        if (debugmsg)
            std::cout<<"Performing "<<ToBeSplitted.size()<< " fold split ops"<<std::endl;

        SplitLev<FaceType> splMd(&ToBeSplitted);
        EdgePred<FaceType> eP(&ToBeSplitted);

        //do the final split

        size_t NumV0=mesh.vert.size();
        bool done=vcg::tri::RefineE<MeshType,SplitLevType,EdgePredType>(mesh,splMd,eP);
        size_t NumV1=mesh.vert.size();

        if (debugmsg)
            std::cout<<"Added "<<NumV1-NumV0<<" vertices"<<std::endl;


        mesh.UpdateDataStructures();
        mesh.SetFeatureFromTable();
        return done;
    }




    static void Perturb(VertexType &v,ScalarType Magnitudo)
    {
        ScalarType eps=std::numeric_limits<ScalarType>::epsilon()*Magnitudo;
        //take a random direction
        size_t granularity=10000;
        int IntX=(rand()%granularity)-granularity/2;
        int IntY=(rand()%granularity)-granularity/2;
        int IntZ=(rand()%granularity)-granularity/2;
        CoordType Dir=CoordType(IntX,IntY,IntZ);
        Dir.Normalize();
        Dir*=eps;
        //std::cout<<Dir.X()<<";"<<Dir.Y()<<";"<<Dir.Z()<<std::endl;
        v.P()+=Dir;
    }

    static size_t NumDuplicatedV(MeshType &mesh)
    {
        std::set<CoordType> Pos;
        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
        size_t numDupl=0;
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (Pos.count(mesh.vert[i].P())>0)
            {
                mesh.vert[i].SetS();
                numDupl++;
            }
            Pos.insert(mesh.vert[i].P());
        }
        return numDupl;
    }

    static bool RepositionDuplicatedV(MeshType &mesh)
    {
        size_t NumD=NumDuplicatedV(mesh);
        if (NumD==0)return false;
        //int dilate_step=0;
        ScalarType Magnitudo=2;
        do
        {
            std::cout<<"Repositioning "<<NumD<<" duplicated vertices"<<std::endl;

            for (size_t i=0;i<mesh.vert.size();i++)
                if (mesh.vert[i].IsS())Perturb(mesh.vert[i],Magnitudo);

            Magnitudo*=2;

            NumD=NumDuplicatedV(mesh);
        }
        while(NumD>0);
        vcg::tri::UpdateBounding<MeshType>::Box(mesh);
        vcg::tri::UpdateNormal<MeshType>::PerVertexNormalizedPerFace(mesh);
        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(mesh);
        return true;
    }

    static bool RemoveZeroAreaF(MeshType &mesh,bool debugmsg=false)
    {
        //        int nonManifV=0;
        //        int degF=0;
        int zeroAFace=0;
        bool modified=false;
        ScalarType Magnitudo=2;
        do{
            modified=false;
            for (size_t i=0;i<mesh.face.size();i++)
            {
                if (vcg::DoubleArea(mesh.face[i])>0)continue;
                Perturb(*mesh.face[i].V(0),Magnitudo);
                Perturb(*mesh.face[i].V(1),Magnitudo);
                Perturb(*mesh.face[i].V(2),Magnitudo);
                modified=true;
                zeroAFace++;
            }
            Magnitudo*=2;
        }while (modified);
        vcg::tri::Allocator<MeshType>::CompactEveryVector(mesh);

        if (debugmsg)
            std::cout<<"Adjusted "<<zeroAFace<<" zero area faces"<<std::endl;
        //        std::cout<<"Removed "<<degF<<" degenerate faces"<<std::endl;
        //        std::cout<<"Removed "<<zeroAFace<<" nonManifV "<<std::endl;
        mesh.UpdateDataStructures();
        return modified;
    }

    static bool RemoveSmallComponents(MeshType &mesh,size_t min_size=10)
    {
        std::vector< std::pair<int, typename MeshType::FacePointer> > CCV;
        vcg::tri::Clean<MeshType>::ConnectedComponents(mesh, CCV);
        vcg::tri::ConnectedComponentIterator<MeshType> ci;
        //vcg::tri::UpdateFlags<MeshType>::FaceClearS(*this);
        bool has_removed=false;
        for(unsigned int i=0;i<CCV.size();++i)
        {
            if (CCV[i].first>min_size)continue;
            has_removed=true;
            //std::vector<typename MeshType::FacePointer> FPV;
            for(ci.start(mesh,CCV[i].second);!ci.completed();++ci)
                vcg::tri::Allocator<MeshType>::DeleteFace(mesh,(*(*ci)));
        }
        if (has_removed)
        {
            vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(mesh);
            vcg::tri::Allocator<MeshType>::CompactEveryVector(mesh);
            mesh.UpdateDataStructures();
        }
        return has_removed;
    }

    static bool SolvePrecisionIssues(MeshType &mesh)
    {
        srand(0);
        bool HasRepositioned=false;
        bool HasSolvedZeroF=false;
        bool HasModified=false;
        do{
            HasRepositioned=RepositionDuplicatedV(mesh);
            HasSolvedZeroF=RemoveZeroAreaF(mesh);
            HasModified|=HasRepositioned;
            HasModified|=HasSolvedZeroF;
        }while (HasRepositioned || HasSolvedZeroF);
        return HasModified;
    }


    //refine the faces with 3 edge selected for field constraints
    static bool RefineInternalFacesStepFromEdgeSel(MeshType &mesh)
    {
        mesh.InitFeatureCoordsTable();
        std::vector<int> to_refine_face;
        for (size_t i=0;i<mesh.face.size();i++)
        {
            //find the number of edges
            int Num=0;
            for (size_t j=0;j<3;j++)
            {
                if (!mesh.face[i].IsFaceEdgeS(j))continue;
                Num++;
            }
            if (Num==3)
                to_refine_face.push_back(i);
        }
        if (to_refine_face.size()==0)return false;

        std::cout<<"Performing "<<to_refine_face.size()<< " face refinement ops"<<std::endl;
        for (size_t j=0;j<to_refine_face.size();j++)
        {
            int IndexF=to_refine_face[j];
            CoordType PD1=mesh.face[IndexF].PD1();
            CoordType PD2=mesh.face[IndexF].PD2();
            CoordType NewPos=(mesh.face[IndexF].P(0)+
                              mesh.face[IndexF].P(1)+
                              mesh.face[IndexF].P(2))/3;
            vcg::tri::Allocator<MeshType>::AddVertex(mesh,NewPos);
            VertexType *V0=mesh.face[IndexF].V(0);
            VertexType *V1=mesh.face[IndexF].V(1);
            VertexType *V2=mesh.face[IndexF].V(2);
            VertexType *V3=&mesh.vert.back();
            mesh.face[IndexF].V(2)=V3;
            vcg::tri::Allocator<MeshType>::AddFace(mesh,V1,V2,V3);
            mesh.face.back().PD1()=PD1;
            mesh.face.back().PD2()=PD2;
            vcg::tri::Allocator<MeshType>::AddFace(mesh,V2,V0,V3);
            mesh.face.back().PD1()=PD1;
            mesh.face.back().PD2()=PD2;
        }
        mesh.UpdateDataStructures();
        mesh.SetFeatureFromTable();
        return true;
    }


    static bool SplitAdjacentEdgeSharpFromEdgeSel(MeshType &mesh)
    {
        typedef SplitLev<FaceType> SplitLevType;
        typedef EdgePred<FaceType> EdgePredType;

        mesh.InitFeatureCoordsTable();
        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
        //InitFaceEdgeSelFromFeatureSeq();

        std::set<std::pair<CoordType,CoordType> > EdgePos;

        for (size_t i=0;i<mesh.face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                if (!mesh.face[i].IsFaceEdgeS(j))continue;
                int VIndex0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                int VIndex1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                CoordType P0=mesh.vert[VIndex0].P();
                CoordType P1=mesh.vert[VIndex1].P();
                mesh.vert[VIndex0].SetS();
                mesh.vert[VIndex1].SetS();
                EdgePos.insert(std::pair<CoordType,CoordType>(std::min(P0,P1),std::max(P0,P1)));
            }
        }

        //then save the edges to be splitted
        std::map<CoordPair,CoordType> ToBeSplitted;
        for (size_t i=0;i<mesh.face.size();i++)
        {
            //find the number of edges
            int Num=0;
            for (size_t j=0;j<3;j++)
            {
                int VIndex0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                int VIndex1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                if ((!mesh.vert[VIndex0].IsS())||(!mesh.vert[VIndex1].IsS()))continue;
                CoordType P0=mesh.vert[VIndex0].P();
                CoordType P1=mesh.vert[VIndex1].P();
                std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
                if (EdgePos.count(Key)==1){Num++;continue;}

                ToBeSplitted[Key]=(P0+P1)/2;
            }
            assert(Num<=2);//this should be already solved
        }
        std::cout<<"Performing "<<ToBeSplitted.size()<< " split ops"<<std::endl;

        SplitLevType splMd(&ToBeSplitted);
        EdgePredType eP(&ToBeSplitted);

        //do the final split
        bool done=vcg::tri::RefineE<MeshType,SplitLevType,EdgePredType >(mesh,splMd,eP);

        mesh.UpdateDataStructures();
        mesh.SetFeatureFromTable();
        return done;
    }



    static bool RemoveFolds(MeshType &mesh,
                            ScalarType MinDot=-0.99,
                            bool debugmsg=false)
    {

        mesh.InitFeatureCoordsTable();

        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
                mesh.face[i].ClearFaceEdgeS(j);

        ScalarType AvgEdge=0;
        size_t Num=0;
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                AvgEdge+=(mesh.face[i].P0(j)-mesh.face[i].P1(j)).Norm();
                Num++;
                if (vcg::face::IsBorder(mesh.face[i],j))continue;
                if (!IsFold(mesh.face[i],j,MinDot))continue;
                mesh.face[i].SetFaceEdgeS(j);
            }
        AvgEdge/=Num;

        size_t NumV0=mesh.vert.size();
        vcg::tri::CutMeshAlongSelectedFaceEdges<MeshType>(mesh);
        mesh.UpdateDataStructures();
        size_t NumV1=mesh.vert.size();

        if (debugmsg)
            std::cout<<"Added "<<NumV1-NumV0<<" vertices"<<std::endl;

        for (size_t i=NumV0;i<NumV1;i++)
            mesh.vert[i].P()+=mesh.vert[i].N()*AvgEdge*0.00001;

        mesh.SetFeatureFromTable();
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<(int)mesh.face[i].VN();j++)
            {
                if (!vcg::face::IsBorder(mesh.face[i],j))continue;
                mesh.face[i].SetFaceEdgeS(j);
            }
        return (NumV1>NumV0);
    }


    static bool RemoveNonManifold(MeshType &mesh)
    {
        bool modified=false;
        ScalarType interval=std::numeric_limits<float>::epsilon()*10;
        int i = 0;
        int removed = 0;
        int max_step=50;
        do
        {
            modified| vcg::tri::Clean<MeshType>::SplitNonManifoldVertex(mesh,interval);
            std::cout << "removed " << removed << " non manifold vertices..." << std::endl;
        }
        while (removed > 0 && i < max_step);

    }

    static bool MakeOrientable(MeshType &mesh)
    {
        bool oriented = false, orientable = false;
        vcg::tri::Clean<MeshType>::OrientCoherentlyMesh(mesh, oriented, orientable);
        mesh.UpdateDataStructures();
        return (!orientable);
    }

    static bool RemoveZeroAreaFaces(MeshType &mesh)
    {
        //REMOVE ZERO AREA FACES
        int cleaned=vcg::tri::Clean<MeshType>::RemoveZeroAreaFace(mesh);
        bool modified=(cleaned>0);

        if (modified)
            mesh.UpdateDataStructures();

        return modified;
    }

    static bool RemoveNonManifolds(MeshType &mesh,
                                   bool debugMsg=false)
    {
        int splitV = 0, openings = 0;
        bool modified=false;
        ScalarType interval=std::numeric_limits<float>::epsilon()*10;
        do
        {
            vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
            splitV = vcg::tri::Clean<MeshType>::SplitNonManifoldVertex(mesh,interval);
            openings = AutoRemesher<MeshType>::openNonManifoldEdges(mesh, interval);
            modified|=(openings>0);

            if (debugMsg)
                std::cout << "Opened " << openings << " non manifold edges" << std::endl;

        } while (splitV > 0 || openings > 0);

        if (modified)
            mesh.UpdateDataStructures();

        return modified;
    }

    static bool RemoveColinearFaces(MeshType & m, const ScalarType colinearThr = 0.001)
    {
        typedef vcg::GridStaticPtr<FaceType, ScalarType> StaticGrid;

        vcg::tri::UpdateTopology<MeshType>::FaceFace(m);

        MeshType projectMesh;
        vcg::tri::Append<MeshType, MeshType>::MeshCopy(projectMesh, m);
        vcg::tri::UpdateBounding<MeshType>::Box(projectMesh);
        vcg::tri::UpdateNormal<MeshType>::PerVertexNormalizedPerFace(projectMesh);
        StaticGrid grid;
        grid.Set(projectMesh.face.begin(), projectMesh.face.end());


        int count = 0;
        int iter = 0;
        bool modified=false;
        do
        {
            vcg::tri::UpdateTopology<MeshType>::FaceFace(m);
            vcg::tri::UnMarkAll(m);

            count = 0;
            for (size_t i = 0; i < size_t(m.FN()); ++i)
            {
                FaceType & f = m.face[i];

                ScalarType quality = vcg::QualityRadii(f.cP(0), f.cP(1), f.cP(2));

                if (quality <= colinearThr)
                {
                    //find longest edge
                    double edges[3];
                    edges[0] = vcg::Distance(f.cP(0), f.cP(1));
                    edges[1] = vcg::Distance(f.cP(1), f.cP(2));
                    edges[2] = vcg::Distance(f.cP(2), f.cP(0));

                    ScalarType smallestEdge = std::min(edges[0], std::min(edges[1], edges[2]));
                    int longestIdx = std::find(edges, edges+3, std::max(std::max(edges[0], edges[1]), edges[2])) - (edges);

                    if (vcg::tri::IsMarked(m, f.V2(longestIdx)))
                        continue;


                    auto f1 = f.cFFp(longestIdx);
                    vcg::tri::Mark(m,f.V2(longestIdx));
                    if (!vcg::face::IsBorder(f, longestIdx) && vcg::face::IsManifold(f, longestIdx) && vcg::face::checkFlipEdgeNotManifold<FaceType>(f, longestIdx))  {

                        // Check if EdgeFlipping improves quality
                        FacePointer g = f.FFp(longestIdx); int k = f.FFi(longestIdx);
                        vcg::Triangle3<ScalarType> t1(f.P(longestIdx), f.P1(longestIdx), f.P2(longestIdx)), t2(g->P(k), g->P1(k), g->P2(k)),
                                t3(f.P(longestIdx), g->P2(k), f.P2(longestIdx)), t4(g->P(k), f.P2(longestIdx), g->P2(k));

                        auto n1 = vcg::TriangleNormal(t1);
                        auto n2 = vcg::TriangleNormal(t2);
                        auto n3 = vcg::TriangleNormal(t3);
                        auto n4 = vcg::TriangleNormal(t4);

                        auto biggestSmallest = vcg::DoubleArea(t1) > vcg::DoubleArea(t2) ? std::make_pair(t1, t2) : std::make_pair(t2, t1);
                        auto areaRatio = vcg::DoubleArea(biggestSmallest.first) / vcg::DoubleArea(biggestSmallest.second);

                        bool normalCheck = true;
                        //                        if (n1.Norm() > 0.001 && n2.Norm() > 0.001)
                        {
                            auto referenceNormal = vcg::NormalizedTriangleNormal(biggestSmallest.first);

                            normalCheck &= vcg::NormalizedTriangleNormal(t3) * referenceNormal >= 0.95;
                            normalCheck &= vcg::NormalizedTriangleNormal(t4) * referenceNormal >= 0.95;
                        }

                        bool areaCheck = false;
                        if (areaRatio > 1000)
                        {
                            areaCheck |= vcg::DoubleArea(t3) / vcg::DoubleArea(biggestSmallest.second) > 1000 && vcg::DoubleArea(t4) / vcg::DoubleArea(biggestSmallest.second) > 1000;
                        }

                        if ((normalCheck) && (areaCheck || std::min( QualityFace(t1), QualityFace(t2) ) <= std::min( QualityFace(t3), QualityFace(t4))))
                        {
                            ScalarType dist;
                            CoordType closest;
                            auto fp0 = vcg::tri::GetClosestFaceBase(projectMesh, grid, vcg::Barycenter(t3), smallestEdge/4., dist, closest);
                            if (fp0 == NULL)
                                continue;

                            auto fp1 = vcg::tri::GetClosestFaceBase(projectMesh, grid, vcg::Barycenter(t4), smallestEdge/4., dist, closest);
                            if (fp1 == NULL)
                                continue;

                            vcg::face::FlipEdgeNotManifold<FaceType>(f, longestIdx);
                            modified=true;
                            ++count;
                        }
                    }
                }
            }
        } while (count && ++iter < 75);
        if (modified)
            m.UpdateDataStructures();

        return modified;
    }

public:

    static bool SolveGeometricArtifactsStep(MeshType &mesh)
    {
        bool modified=false;

//        //REMOVE COLLINEAR FACES
        modified|=RemoveColinearFaces(mesh);

        //REMOVE ZERO AREA FACES
        modified|=RemoveZeroAreaFaces(mesh);

        //SPLIT NON MANIFOLD FACES
        modified|=RemoveNonManifolds(mesh);

        //REMOVE MINIMAL CONNECTED COMPONENTS
        modified|=RemoveSmallComponents(mesh);

        //MAKE ORIENTABLE
        modified|=MakeOrientable(mesh);

        //SOLVE POSSIBLE PRECISION ISSUES
        modified|=SolvePrecisionIssues(mesh);

//        //THEN SPLIT OR REMOVE 180 FOLDS
        modified|=SplitFolds(mesh);
        modified|=RemoveFolds(mesh);

//        //REMOVE POSSIBLE PRECISION ISSUES
        modified|=SolvePrecisionIssues(mesh);

        return modified;
    }

public:

    static void SolveGeometricArtifacts(MeshType &mesh,size_t max_steps=10)
    {
        size_t currS=0;
        while ((SolveGeometricArtifactsStep(mesh))&&(currS<max_steps))
            currS++;

        //AutoRemesher<MeshType>::collapseSurvivingMicroEdges(mesh,0.001, const ScalarType edgeRatio = 0.025);

        vcg::tri::Allocator<MeshType>::CompactEveryVector(mesh);
        mesh.UpdateDataStructures();
    }

    static void RefineIfNeeded(MeshType &mesh)
    {

        mesh.UpdateDataStructures();
        bool has_refined=false;
        do
        {
            has_refined=false;
            has_refined|=RefineInternalFacesStepFromEdgeSel(mesh);
            //SolvePrecisionIssues(mesh);//CHECK
            has_refined|=SplitAdjacentEdgeSharpFromEdgeSel(mesh);
            //SolvePrecisionIssues(mesh);//CHECK

            //has_refined|=SplitAdjacentTerminalVertices();
            //has_refined|=SplitEdgeSharpSharingVerticesFromEdgeSel();
        }while (has_refined);
        mesh.InitEdgeType();
    }

    static void InitSharpFeatures(MeshType &mesh,
                                  ScalarType sharp_feature_thr,
                                  size_t feature_erode_dilate)
    {
        mesh.UpdateDataStructures();
        mesh.InitSharpFeatures(sharp_feature_thr);
        mesh.ErodeDilate(feature_erode_dilate);
    }

    struct BatchParam
    {
        bool UpdateSharp=true;
        bool DoRemesh=true;
        bool surf_dist_check=true;
        ScalarType sharp_feature_thr=35;
        size_t feature_erode_dilate=4;
        ScalarType SharpFactor=6;
        size_t remesher_iterations=15;
        ScalarType remesher_aspect_ratio=0.3;
        ScalarType remesher_termination_delta = 10000;
    };

    static void BatchProcess(MeshType &mesh,BatchParam &BPar,
                             typename vcg::tri::FieldSmoother<MeshType>::SmoothParam &FieldParam)
    {
        mesh.UpdateDataStructures();

        // SELECT SHARP FEATURES
        if (BPar.UpdateSharp)
            MeshPrepocess<MeshType>::InitSharpFeatures(mesh,BPar.sharp_feature_thr,BPar.feature_erode_dilate);

//        MeshPrepocess<MeshType>::SolveGeometricArtifacts(mesh);

//        // SELECT SHARP FEATURES
//        if (BPar.UpdateSharp)
//            MeshPrepocess<MeshType>::InitSharpFeatures(mesh,BPar.sharp_feature_thr,BPar.feature_erode_dilate);

        //DO REMESH IF NEEDED
        if (BPar.DoRemesh)
        {
            typename AutoRemesher<MeshType>::Params RemPar;
            RemPar.iterations   = BPar.remesher_iterations;
            RemPar.targetAspect = BPar.remesher_aspect_ratio;
            RemPar.targetDeltaFN= BPar.remesher_termination_delta;
            RemPar.surfDistCheck = BPar.surf_dist_check;

            //AutoRemesher<MeshType>::Remesh2(mesh,RemPar);
            AutoRemesher<MeshType>::RemeshAdapt(mesh,RemPar);

        }
        mesh.InitFeatureCoordsTable();

//        //THEN UPDATE SHARP FEATURES
//        if (BPar.UpdateSharp)
//            MeshPrepocess<MeshType>::InitSharpFeatures(mesh,BPar.sharp_feature_thr,BPar.feature_erode_dilate);


        //SOLVE POSSIBLE GEOMETRIC ARTIFACTS AFTER REFINEMENT
        MeshPrepocess<MeshType>::SolveGeometricArtifacts(mesh);

        //REFINE THE MESH IF NEEDED TO BE CONSISTENT WHEN COMPUTING FIELD
        MeshPrepocess<MeshType>::RefineIfNeeded(mesh);

//        //THEN UPDATE SHARP FEATURES
//        if (BPar.UpdateSharp)
//            MeshPrepocess<MeshType>::InitSharpFeatures(mesh,BPar.sharp_feature_thr,BPar.feature_erode_dilate);
//        mesh.ErodeDilate(BPar.feature_erode_dilate);

        MeshFieldSmoother<MeshType>::AutoSetupParam(mesh,FieldParam,BPar.SharpFactor);

        //THEN SMOOTH THE FIELD
        std::cout << "[fieldComputation] Smooth Field Computation..." << std::endl;
        MeshFieldSmoother<MeshType>::SmoothField(mesh,FieldParam);
    }

    static void SaveAllData(MeshType &tri_mesh,const std::string &pathM)
    {
        std::string projM=pathM;
        size_t indexExt=projM.find_last_of(".");
        projM=projM.substr(0,indexExt);
        std::string meshName=projM+std::string("_rem.obj");
        std::string fieldName=projM+std::string("_rem.rosy");
        std::string sharpName=projM+std::string("_rem.sharp");
        std::cout<<"Saving Mesh TO:"<<meshName.c_str()<<std::endl;
        std::cout<<"Saving Field TO:"<<fieldName.c_str()<<std::endl;
        std::cout<<"Saving Sharp TO:"<<sharpName.c_str()<<std::endl;
        tri_mesh.SaveTriMesh(meshName.c_str());
        tri_mesh.SaveField(fieldName.c_str());
        tri_mesh.SaveSharpFeatures(sharpName.c_str());
    }

};




#endif
