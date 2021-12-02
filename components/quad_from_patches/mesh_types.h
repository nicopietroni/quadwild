#ifndef DEFAULTMESHTYPES_H
#define DEFAULTMESHTYPES_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/clean.h>
#include <wrap/io_trimesh/export.h>

/* ----- Polygon mesh ----- */

class PolyVertex;
class PolyFace;
class PolyEdge;

struct MyPolyTypes : public vcg::UsedTypes<
        vcg::Use<PolyVertex>::AsVertexType,
        vcg::Use<PolyEdge>::AsEdgeType,
        vcg::Use<PolyFace>::AsFaceType>{};

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

class PolyVertex : public vcg::Vertex<MyPolyTypes,
        vcg::vertex::Coord3f,
        vcg::vertex::Normal3f,
        vcg::vertex::Color4b,
        vcg::vertex::Qualityf,
        vcg::vertex::BitFlags,
        vcg::vertex::VFAdj,
        vcg::vertex::CurvatureDirf>{};

class PolyFace : public vcg::Face<
        MyPolyTypes,
        vcg::face::PolyInfo,
        vcg::face::VertexRef,
        vcg::face::Normal3f,
        vcg::face::Color4b,
        vcg::face::Qualityf,
        vcg::face::BitFlags,
        vcg::face::PFVAdj,
        vcg::face::PFFAdj,
        vcg::face::PVFAdj,
        vcg::face::CurvatureDirf,
        vcg::face::Mark> {};

class PolyEdge : public vcg::Edge<
        MyPolyTypes,
        vcg::edge::VertexRef,
        vcg::edge::BitFlags> {};

class PolyMesh : public vcg::tri::TriMesh<
        std::vector<PolyVertex>,
        std::vector<PolyEdge>,
        std::vector<PolyFace>> {};


/* ----- Triangle mesh ----- */

class TriangleVertex;
class TriangleFace;
struct MyTriangleTypes : public vcg::UsedTypes<
        vcg::Use<TriangleVertex>::AsVertexType,
        vcg::Use<TriangleFace>::AsFaceType>{};

class TriangleVertex : public vcg::Vertex<
        MyTriangleTypes,
        vcg::vertex::Coord3f,
        vcg::vertex::Normal3f,
        vcg::vertex::VFAdj,
        vcg::vertex::Color4b,
        vcg::vertex::Qualityf,
        vcg::vertex::BitFlags,
        vcg::vertex::CurvatureDirf>{};

class TriangleFace : public vcg::Face<
        MyTriangleTypes,
        vcg::face::VertexRef,
        vcg::face::Normal3f,
        vcg::face::Color4b,
        vcg::face::Qualityd,
        vcg::face::BitFlags,
        vcg::face::FFAdj,
        vcg::face::VFAdj,
        vcg::face::CurvatureDirf,
        vcg::face::Mark> {};

class TriangleMesh : public vcg::tri::TriMesh<
        std::vector<TriangleVertex>,
        std::vector<TriangleFace> > {



//    void Perturb(VertexType &v,ScalarType Magnitudo)
//    {
//        ScalarType eps=std::numeric_limits<ScalarType>::epsilon()*Magnitudo;
//        //take a random direction
//        size_t granularity=10000;
//        int IntX=(rand()%granularity)-granularity/2;
//        int IntY=(rand()%granularity)-granularity/2;
//        int IntZ=(rand()%granularity)-granularity/2;
//        CoordType Dir=CoordType(IntX,IntY,IntZ);
//        Dir.Normalize();
//        Dir*=eps;
//        //std::cout<<Dir.X()<<";"<<Dir.Y()<<";"<<Dir.Z()<<std::endl;
//        v.P()+=Dir;
//    }

//    size_t NumDuplicatedV()
//    {
//        std::set<CoordType> Pos;
//        vcg::tri::UpdateSelection<MeshType>::VertexClear(*this);
//        size_t numDupl=0;
//        for (size_t i=0;i<vert.size();i++)
//        {
//            if (Pos.count(vert[i].P())>0)
//            {
//                vert[i].SetS();
//                numDupl++;
//            }
//            Pos.insert(vert[i].P());
//        }
//        return numDupl;
//    }

//    bool RepositionDuplicatedV()
//    {
//        size_t NumD=NumDuplicatedV();
//        std::cout<<"There are "<<NumD<<" duplicated vertices"<<std::endl;

//        if (NumD==0)return false;
//        //int dilate_step=0;
//        ScalarType Magnitudo=2;
//        do
//        {
//            std::cout<<"Repositioning "<<NumD<<" duplicated vertices"<<std::endl;


//            //            dilate_step++;
//            //            for (size_t i=0;i<dilate_step;i++)
//            //            {
//            //                vcg::tri::UpdateSelection<MeshType>::FaceFromVertexLoose(*this);
//            //                vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(*this);
//            //            }
//            for (size_t i=0;i<vert.size();i++)
//                if (vert[i].IsS())Perturb(vert[i],Magnitudo);

//            Magnitudo*=2;
//            //vcg::tri::Smooth<MeshType>::VertexCoordLaplacian(*this,1,true);
//            NumD=NumDuplicatedV();
//        }
//        while(NumD>0);
//        vcg::tri::UpdateBounding<MeshType>::Box(*this);
//        vcg::tri::UpdateNormal<MeshType>::PerVertexNormalizedPerFace(*this);
//        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(*this);
//        return true;
//    }

//    bool RemoveZeroAreaF()
//    {
//        //        int nonManifV=0;
//        //        int degF=0;

//        int zeroAFace=0;
//        bool modified=false;
//        ScalarType Magnitudo=2;
//        do{
//            modified=false;
//            for (size_t i=0;i<face.size();i++)
//            {
//                if (vcg::DoubleArea(face[i])>0)continue;
//                Perturb(*face[i].V(0),Magnitudo);
//                Perturb(*face[i].V(1),Magnitudo);
//                Perturb(*face[i].V(2),Magnitudo);
//                modified=true;
//                zeroAFace++;
//            }
//            Magnitudo*=2;
//        }while (modified);
//        vcg::tri::Allocator<MeshType>::CompactEveryVector(*this);
//        std::cout<<"Adjusted "<<zeroAFace<<" zero area faces"<<std::endl;
//        //        std::cout<<"Removed "<<degF<<" degenerate faces"<<std::endl;
//        //        std::cout<<"Removed "<<zeroAFace<<" nonManifV "<<std::endl;
//        UpdateAttributes();
//        return modified;
//    }


public:

    void UpdateAttributes()
    {
        vcg::tri::UpdateNormal<TriangleMesh>::PerFaceNormalized(*this);
        vcg::tri::UpdateNormal<TriangleMesh>::PerVertexNormalized(*this);
        vcg::tri::UpdateBounding<TriangleMesh>::Box(*this);
        vcg::tri::UpdateTopology<TriangleMesh>::FaceFace(*this);
        vcg::tri::UpdateTopology<TriangleMesh>::VertexFace(*this);
        vcg::tri::UpdateFlags<TriangleMesh>::FaceBorderFromFF(*this);
        vcg::tri::UpdateFlags<TriangleMesh>::VertexBorderFromFaceBorder(*this);
    }

//    void SolveGeometricIssues()
//    {
//        srand(0);
//        bool HasRepositioned=false;
//        bool HasSolvedZeroF=false;
//        do{
//            HasRepositioned=RepositionDuplicatedV();
//            HasSolvedZeroF=RemoveZeroAreaF();
//        }while (HasRepositioned || HasSolvedZeroF);
//        UpdateAttributes();

//    }
};

//void SelectMeshPatchBorders(TriangleMesh &mesh,const std::vector<size_t>  &FacePatches)
//{
//    //std::set<std::pair<size_t,size_t> > BorderPatches;
//    vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);

//    assert(FacePatches.size()==mesh.face.size());
//    //first add borders
//    for (size_t i=0;i<mesh.face.size();i++)
//        for (size_t j=0;j<3;j++)
//        {
//            if (mesh.face[i].IsB(j))
//            {
//                    mesh.face[i].SetFaceEdgeS(j);
//            }
//            else
//            {
//                size_t FOpp=vcg::tri::Index(mesh,mesh.face[i].FFp(j));
//                assert(FOpp!=i);
//                if (FacePatches[i]<0)assert(FacePatches[i]==-1);
//                //assert(FacePatches[i]>=0);

//                if (FacePatches[i]!=FacePatches[FOpp])
//                {
//                        mesh.face[i].SetFaceEdgeS(j);
//                }
//            }
//        }
//}

void ComputePerFacePatch(TriangleMesh &mesh,
                         const std::vector<std::vector<size_t> >  &PatchFaces,
                         std::vector<size_t>  &FacePatches)
{
    FacePatches.clear();
    FacePatches.resize(mesh.face.size(),mesh.face.size()+1);
    for (size_t i=0;i<PatchFaces.size();i++)
        for (size_t j=0;j<PatchFaces[i].size();j++)
        {
            size_t IndexF=PatchFaces[i][j];
            assert(IndexF<mesh.face.size());
            FacePatches[IndexF]=i;
        }
    for (size_t i=0;i<FacePatches.size();i++)
        assert(FacePatches[i]<mesh.face.size());
}

void SelectVertMeshPatchBorders(TriangleMesh &mesh,const std::vector<size_t>  &FacePatches)
{
    //std::set<std::pair<size_t,size_t> > BorderPatches;
    vcg::tri::UpdateSelection<TriangleMesh>::VertexClear(mesh);

    assert(FacePatches.size()==mesh.face.size());
    //first add borders
    for (size_t i=0;i<mesh.face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            if (mesh.face[i].IsB(j))
            {
                mesh.face[i].V0(j)->SetS();
                mesh.face[i].V1(j)->SetS();
            }
            else
            {
                size_t FOpp=vcg::tri::Index(mesh,mesh.face[i].FFp(j));
                assert(FOpp!=i);
                //if (FacePatches[i]<0)assert(FacePatches[i]==-1);
                //assert(FacePatches[i]>=0);

                if (FacePatches[i]!=FacePatches[FOpp])
                {
                    mesh.face[i].V0(j)->SetS();
                    mesh.face[i].V1(j)->SetS();
                }
            }
        }
}

void CheckVertInPatch(TriangleMesh &mesh,
                      const std::vector<size_t>  &PatchFaces,
                      const std::vector<size_t>  &PatchCorners)
{
    //std::set<std::pair<size_t,size_t> > BorderPatches;
    vcg::tri::UpdateSelection<TriangleMesh>::VertexClear(mesh);

    //select patch vertices
    for (size_t i=0;i<PatchFaces.size();i++)
    {
        size_t IndexF=PatchFaces[i];
        for (size_t j=0;j<3;j++)
            mesh.face[IndexF].V(j)->SetS();
    }

    for (size_t i=0;i<PatchCorners.size();i++)
    {
        size_t IndexV=PatchCorners[i];
        assert(mesh.vert[IndexV].IsS());
    }
}

void CheckManifoldPatch(TriangleMesh &mesh,
                      const std::vector<size_t>  &PatchFaces)
{
    TriangleMesh testM;
    vcg::tri::UpdateSelection<TriangleMesh>::Clear(mesh);

    //select patch vertices
    for (size_t i=0;i<PatchFaces.size();i++)
    {
        size_t IndexF=PatchFaces[i];
        mesh.face[IndexF].SetS();
    }
    vcg::tri::UpdateSelection<TriangleMesh>::VertexFromFaceLoose(mesh);
    vcg::tri::Append<TriangleMesh,TriangleMesh>::Mesh(testM,mesh,true);
    testM.UpdateAttributes();
    assert(vcg::tri::Clean<TriangleMesh>::CountNonManifoldEdgeFF(testM)==0);
    assert(vcg::tri::Clean<TriangleMesh>::CountNonManifoldVertexFF(testM)==0);
}

void ExportPatch(TriangleMesh &mesh,
                 const std::vector<size_t>  &PatchFaces)
{
    TriangleMesh testM;
    vcg::tri::UpdateSelection<TriangleMesh>::Clear(mesh);

    //select patch vertices
    for (size_t i=0;i<PatchFaces.size();i++)
    {
        size_t IndexF=PatchFaces[i];
        mesh.face[IndexF].SetS();
    }
    vcg::tri::UpdateSelection<TriangleMesh>::VertexFromFaceLoose(mesh);
    vcg::tri::Append<TriangleMesh,TriangleMesh>::Mesh(testM,mesh,true);
    testM.UpdateAttributes();
    assert(vcg::tri::Clean<TriangleMesh>::CountNonManifoldEdgeFF(testM)==0);
    assert(vcg::tri::Clean<TriangleMesh>::CountNonManifoldVertexFF(testM)==0);
    vcg::tri::io::ExporterPLY<TriangleMesh>::Save(testM,"test.ply");
}

void CheckIntegrity(TriangleMesh &mesh,
                    const std::vector<std::vector<size_t> >  &PatchFaces,
                    const std::vector<std::vector<size_t>> &PatchCorners)
{
    std::vector<size_t>  FacePatches;
    ComputePerFacePatch(mesh,PatchFaces,FacePatches);
    SelectVertMeshPatchBorders(mesh,FacePatches);
    for (size_t i=0;i<PatchCorners.size();i++)
        for (size_t j=0;j<PatchCorners[i].size();j++)
        {
            size_t IndexV=PatchCorners[i][j];
            assert(mesh.vert[IndexV].IsS());
        }
    //then check the corners belong to the patch
    for (size_t i=0;i<PatchFaces.size();i++)
    {
//        if (i==41)
//        {
//            ExportPatch(mesh,PatchFaces[i]);
//        }
        //std::cout<<"Checking Partition: "<<i<<std::endl;
        //CheckManifoldPatch(mesh,PatchFaces[i]);
        CheckVertInPatch(mesh,PatchFaces[i],PatchCorners[i]);
    }

}

void OrientIfNeeded(TriangleMesh &mesh,
                    std::vector<std::vector<size_t> > &trimeshPartitions,
                    std::vector<std::vector<size_t>> &trimeshCorners,
                    std::vector<std::pair<size_t,size_t> > &trimeshFeatures,
                    std::vector<size_t> &trimeshFeaturesC)
{
    std::cout<<"CHECKING ORIENTATION"<<std::endl;

    mesh.UpdateAttributes();
    bool IsOriented,IsOrientable;
    IsOriented=vcg::tri::Clean<TriangleMesh>::IsCoherentlyOrientedMesh(mesh);
    if (IsOriented)return;
    std::cout<<"MESH NEED TO BE RE-ORIENTED"<<std::endl;

    //save the index of original v and faces
    for (size_t i=0;i<mesh.vert.size();i++)
        mesh.vert[i].Q()=i;
    for (size_t i=0;i<mesh.face.size();i++)
        mesh.face[i].Q()=i;

    //save sharp features as pair of vertices
    std::set<std::pair<size_t,size_t> > SharpV;
    for (size_t i=0;i<trimeshFeatures.size();i++)
    {
        size_t indexF=trimeshFeatures[i].first;
        size_t indexE=trimeshFeatures[i].second;
        assert(indexE<3);
        assert(indexF<mesh.face.size());
        size_t IndexV0=vcg::tri::Index(mesh,mesh.face[indexF].V0(indexE));
        size_t IndexV1=vcg::tri::Index(mesh,mesh.face[indexF].V1(indexE));
        std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
        SharpV.insert(key);
    }
    //then orient it coherently
    vcg::tri::Clean<TriangleMesh>::OrientCoherentlyMesh(mesh, IsOriented,IsOrientable);
    mesh.UpdateAttributes();

    //then update sharp features
    trimeshFeatures.clear();
    for (size_t i=0;i<mesh.face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            size_t IndexV0=mesh.face[i].V0(j)->Q();
            size_t IndexV1=mesh.face[i].V1(j)->Q();
            assert(IndexV0!=IndexV1);
            assert(IndexV0<mesh.vert.size());
            assert(IndexV1<mesh.vert.size());
            std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            if (SharpV.count(key)>0)
                trimeshFeatures.push_back(std::pair<size_t,size_t>(i,j));
        }

    if (!IsOrientable)
    {
        std::cout<<"MESH IS NOT ORIENTABLE "<<std::endl;
        assert(0);
    }
}

#endif // DEFAULTMESHTYPES_H
