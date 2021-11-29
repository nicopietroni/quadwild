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

#ifndef LOCAL_PARAM_SMOOTH
#define LOCAL_PARAM_SMOOTH

#include <vcg/space/triangle2.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/parametrization/uv_utils.h>
#include <wrap/igl/lscm_parametrization.h>
#include <wrap/io_trimesh/export.h>
#include <vcg/space/distance3.h>


template <class MeshType>
class Local_Param_Smooth
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename vcg::face::Pos<FaceType> PosType;
    typedef typename vcg::Point2<ScalarType> Point2Type;
    typedef typename vcg::Triangle2<ScalarType> Triangle2Type;

public:

    struct UVSmoothParam
    {
        bool FixSel;
        bool LineEdgeSel;
        ScalarType Damp;
        size_t Steps;

        UVSmoothParam()
        {
            FixSel=false;
            LineEdgeSel=true;
            Damp=0.2;
            Steps=5;
        }
    };

private:


    static PosType GetStarInternalPos(MeshType &sub_mesh)
    {
        //get the first pos to start iterating
        int indexE=-1;
        if (!sub_mesh.face[0].V(0)->IsB())indexE=0;
        if (!sub_mesh.face[0].V(1)->IsB())indexE=1;
        if (!sub_mesh.face[0].V(2)->IsB())indexE=2;
        assert(indexE>=0);
        PosType VertPos(&sub_mesh.face[0],indexE);
        return VertPos;
    }

    static void SetStarInternalPos(MeshType &sub_mesh,const CoordType &Pos)
    {
        //get the first pos to start iterating
        int indexE=-1;
        if (!sub_mesh.face[0].V(0)->IsB())indexE=0;
        if (!sub_mesh.face[0].V(1)->IsB())indexE=1;
        if (!sub_mesh.face[0].V(2)->IsB())indexE=2;
        assert(indexE>=0);
        sub_mesh.face[0].P(indexE)=Pos;
    }

    static bool InterpolateUV(MeshType &mesh,
                              const vcg::Point2<ScalarType> &UV,
                              size_t &IndexF,CoordType &bary)
    {
        for (size_t i=0;i<mesh.face.size();i++)
        {
            vcg::Point2<ScalarType> UV0,UV1,UV2;
            UV0=mesh.face[i].V(0)->T().P();
            UV1=mesh.face[i].V(1)->T().P();
            UV2=mesh.face[i].V(2)->T().P();
            bool inside=InterpolationParameters2(UV0,UV1,UV2,UV,bary);
            if (inside)
            {
                IndexF=i;
                return true;
            }
        }
        return false;
    }


    static Point2Type ComputeCentralUV(std::vector<PosType> &StarFPos)
    {
        ScalarType WSum=0;
        Point2Type  CenterUV=Point2Type(0,0);
        for (size_t i=0;i<StarFPos.size();i++)
        {
            CoordType P1=StarFPos[i].VFlip()->P();
            Point2Type UV1=StarFPos[i].VFlip()->T().P();

            size_t next_i=(i+1)%StarFPos.size();

            CoordType P2=StarFPos[next_i].VFlip()->P();
            Point2Type UV2=StarFPos[next_i].VFlip()->T().P();

            CoordType P0=StarFPos[i].V()->P();

            ScalarType angle1=vcg::Angle(P2 - P1, P0 - P1);
            ScalarType angle2=vcg::Angle(P1 - P2, P0 - P2);
            ScalarType weight1 = tan((M_PI * 0.5) - angle1);
            ScalarType weight2 = tan((M_PI * 0.5) - angle2);

            CenterUV+=UV1*weight2;
            CenterUV+=UV2*weight1;

            WSum+=weight1;
            WSum+=weight2;
        }
        CenterUV/=WSum;
        return CenterUV;
    }

    static void ParametrizeCentralUV(MeshType &StarMesh)
    {
        PosType StartVertPos=GetStarInternalPos(StarMesh);
        //get the star around
        std::vector<PosType> StarFPos;
        vcg::face::VFOrderedStarFF(StartVertPos,StarFPos);
        assert(!StartVertPos.V()->IsB());
        assert(!StartVertPos.V()->IsS());
        StartVertPos.V()->T().P()=ComputeCentralUV(StarFPos);
    }

    static bool FindSmoothedPos(MeshType &StarMesh,
                                ScalarType Damp,
                                CoordType &InterpPos)
    {

        ParametrizeCentralUV(StarMesh);

        PosType InternalPos=GetStarInternalPos(StarMesh);
        //vcg::Point2<ScalarType> UVCenter=InternalPos.V()->T().P();

        vcg::Point2<ScalarType> CenterUVTarget(0,0);
        vcg::Point2<ScalarType> SamplePos=CenterUVTarget;//UVCenter*Damp+CenterUVTarget*(1-Damp);

        size_t IndexF;
        CoordType bary;
        bool Interpolated=InterpolateUV(StarMesh,SamplePos,IndexF,bary);
        if (!Interpolated)
        {
            //std::cout<<"No Interp"<<std::endl;
            //bool NeedC;
            ParametrizeBorders(StarMesh);//,false,NeedC);//NeedConformal[i]);
            ParametrizeCentralUV(StarMesh);
            Interpolated=InterpolateUV(StarMesh,SamplePos,IndexF,bary);
            if (!Interpolated)
            {
                vcg::tri::io::ExporterPLY<MeshType>::Save(StarMesh,"test.ply");
                assert(0);
            }
            //return false;
        }

        CoordType P0=StarMesh.face[IndexF].P(0);
        CoordType P1=StarMesh.face[IndexF].P(1);
        CoordType P2=StarMesh.face[IndexF].P(2);
        InterpPos=CoordType(P0*bary.X()+P1*bary.Y()+P2*bary.Z());
        return true;

    }

    static CoordType ProjectOnBasis(const CoordType &testP,
                                    const std::vector<PosType> &ProjBasis)
    {
        if (ProjBasis.size()==0)return testP;
        CoordType closestP=testP;
        ScalarType currD=std::numeric_limits<ScalarType>::max();
        //std::cout<<"ProjBasis Size"<<ProjBasis.size();
        for (size_t i=0;i<ProjBasis.size();i++)
        {
            CoordType P0=ProjBasis[i].V()->P();
            CoordType P1=ProjBasis[i].VFlip()->P();
            vcg::Segment3<ScalarType> S3(P0,P1);
            CoordType clos;
            ScalarType distT;
            vcg::SegmentPointDistance(S3,testP,clos,distT);
            if (distT>currD)continue;
            currD=distT;
            closestP=clos;
        }
        return closestP;
    }

    static void SmoothStep(MeshType &mesh,
                           const std::vector<MeshType*> &StarMeshes,
                           const std::vector<std::vector<int> > &OrigIndex,
                           const std::vector<std::vector<PosType> > &ProjBasis,
                           const UVSmoothParam &UVP,
                           const std::vector<bool> &IsS)
    {
        std::vector<CoordType> TargetP(mesh.vert.size(),CoordType(0,0,0));
        std::vector<size_t> TargetN(mesh.vert.size(),0);

        assert(IsS.size()==mesh.vert.size());
        assert(StarMeshes.size()==mesh.vert.size());
        //assert(NeedConformal.size()==StarMeshes.size());

        //compute target pos just to use for borders
        for (size_t i=0;i<mesh.face.size();i++)
        {
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                VertexType *v0=mesh.face[i].V(j);
                VertexType *v1=mesh.face[i].V((j+1)%mesh.face[i].VN());
                //check to sum up only once
                if ((!vcg::face::IsBorder(mesh.face[i],j))&&(v1<v0))continue;
                size_t IndexV0=vcg::tri::Index(mesh,v0);
                size_t IndexV1=vcg::tri::Index(mesh,v1);
                if (v0->IsB())
                {
                    TargetP[IndexV0]+=v1->P();
                    TargetN[IndexV0]++;
                }
                if (v1->IsB())
                {
                    TargetP[IndexV1]+=v0->P();
                    TargetN[IndexV1]++;
                }
            }
        }

        for (size_t i=0;i<TargetP.size();i++)
        {
            if (!mesh.vert[i].IsB())continue;
            assert(TargetN[i]>0);
            TargetP[i]/=TargetN[i];
        }

        for (size_t i=0;i<StarMeshes.size();i++)
        {
            //if (StarMeshes[i]==NULL)
            if (mesh.vert[i].IsB())
            {
                //assert(StarMeshes[i]==NULL);
                continue;
            }
            else
            {
                assert(StarMeshes[i]!=NULL);
                //compute the pre-parametrization
                //bool NeedC=true;
                ParametrizeBorders((*StarMeshes[i]));//,UVP.LineEdgeSel,NeedC);//NeedConformal[i]);

                CoordType InterpPos;
                //bool Interpolated=FindSmoothedPos((*StarMeshes[i]),NeedC,UVP.Damp,InterpPos);
                bool Interpolated=FindSmoothedPos((*StarMeshes[i]),UVP.Damp,InterpPos);
                if (Interpolated)
                    TargetP[i]=InterpPos;
                else
                    TargetP[i]=mesh.vert[i].P();
            }
        }

        //update mesh position
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (UVP.FixSel && IsS[i])continue;

            mesh.vert[i].P()=mesh.vert[i].P()*UVP.Damp+TargetP[i]*(1-UVP.Damp);

            if (UVP.LineEdgeSel)
                mesh.vert[i].P()=ProjectOnBasis(mesh.vert[i].P(),ProjBasis[i]);
        }

        //update submesh positions
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (StarMeshes[i]==NULL)continue;
            assert(OrigIndex[i].size()==(*StarMeshes[i]).vert.size());
            for (size_t j=0;j<(*StarMeshes[i]).vert.size();j++)
            {
                size_t IndexO=OrigIndex[i][j];
                (*StarMeshes[i]).vert[j].P()=mesh.vert[IndexO].P();
            }
        }

    }



    static void DeriveRegularStarPos(int N,std::vector<Point2Type> &StarPos)
    {
        std::vector<CoordType> TemplatePos;
        vcg::getBaseTemplatePolygon(N,TemplatePos);
        for (size_t i=0;i<TemplatePos.size();i++)
            StarPos.push_back(Point2Type(TemplatePos[i].X(),TemplatePos[i].Y()));
    }

    static void ParametrizeBorders(MeshType &sub_mesh)/*,
                                   const bool UseSharp),
                                   bool &NeedConformal)*/
    {
        PosType StartVertPos=GetStarInternalPos(sub_mesh);

        //get the star around
        std::vector<PosType> StarFPos;
        vcg::face::VFOrderedStarFF(StartVertPos,StarFPos);

        //then get the indexes of the one having sharp edges
        std::vector<size_t> SelIndexes;

        for (size_t i=0;i<StarFPos.size();i++)
            SelIndexes.push_back(i);

        //then simply parametrize the ones that should be fixed on border
        std::vector<Point2Type> UV_Boundary;
        DeriveRegularStarPos(SelIndexes.size(),UV_Boundary);
        vcg::tri::UpdateSelection<MeshType>::VertexClear(sub_mesh);

        for (size_t i=0;i<SelIndexes.size();i++)
        {
            //get the boundary fixed indexes
            size_t currI=SelIndexes[i];
            //set the coordinates
            StarFPos[currI].VFlip()->T().P()=UV_Boundary[i];
            //then select it as it shoould be preserved during param
            StarFPos[currI].VFlip()->SetS();
        }
    }

    static void GetStarSubMeshes(MeshType &mesh,
                                 const PosType &StartVertPos,
                                 MeshType &sub_mesh,
                                 std::vector<PosType> &ProjBasis,
                                 std::vector<int> &OrigIndex)
    {
        OrigIndex.clear();

        std::vector<PosType> StarFPos;
        vcg::face::VFOrderedStarFF(StartVertPos,StarFPos);

        vcg::tri::UpdateSelection<MeshType>::FaceClear(mesh);
        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);

        for (size_t k=0;k<StarFPos.size();k++)
            StarFPos[k].F()->SetS();

        sub_mesh.Clear();
        vcg::tri::Append<MeshType,MeshType>::Mesh(sub_mesh,mesh,true);

        //copy the index
        for (size_t i=0;i<sub_mesh.vert.size();i++)
            OrigIndex.push_back(sub_mesh.vert[i].Q());

        vcg::tri::UpdateTopology<MeshType>::FaceFace(sub_mesh);
        vcg::tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(sub_mesh);

        //then retrieve a pos from the central vert
        bool found=false;
        PosType StartPos1;
        size_t IndexCentral=vcg::tri::Index(mesh,StartVertPos.V());
        for (size_t i=0;i<sub_mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                size_t IndexV=vcg::tri::Index(sub_mesh,sub_mesh.face[i].V(j));
                if (OrigIndex[IndexV]==IndexCentral)
                {
                    found=true;
                    StartPos1=PosType(&sub_mesh.face[i],j);
                    break;
                }
            }
        std::vector<PosType> Star1FPos;
        vcg::face::VFOrderedStarFF(StartPos1,Star1FPos);

        for (size_t k=0;k<Star1FPos.size();k++)
        {
            if (Star1FPos[k].IsBorder() || Star1FPos[k].IsEdgeS())
                ProjBasis.push_back(Star1FPos[k]);
        }

    }

    static void GetStarMeshes(MeshType &mesh,
                              std::vector<MeshType*> &StarMeshes,
                              std::vector<std::vector<PosType> > &ProjBasis,
                              std::vector<std::vector<int> > &OrigIndex)
    {
        //first save the idx of original vertx on quality
        for (size_t i=0;i<mesh.vert.size();i++)
            mesh.vert[i].Q()=i;

        //then create the submeshes
        StarMeshes=std::vector<MeshType*>(mesh.vert.size(),NULL);
        OrigIndex.clear();
        OrigIndex.resize(mesh.vert.size());

        ProjBasis.clear();
        ProjBasis.resize(mesh.vert.size());

        //NeedConformal=std::vector<bool>(mesh.vert.size(),true);
        vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (mesh.face[i].V(j)->IsV())continue;
                //if (mesh.face[i].V(j)->IsB())continue;
                mesh.face[i].V(j)->SetV();

                //allocate the submesh
                size_t IndexV=vcg::tri::Index(mesh,mesh.face[i].V(j));
                assert(StarMeshes[IndexV]==NULL);
                StarMeshes[IndexV]=new MeshType();

                //get the starting pos
                PosType VertPos(&mesh.face[i],j);

                //copy the submesh
                GetStarSubMeshes(mesh,VertPos,(*StarMeshes[IndexV]),ProjBasis[IndexV],OrigIndex[IndexV]);

            }
        vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
    }

public:

    static void Smooth(MeshType &mesh,UVSmoothParam UVP=UVSmoothParam())
    {
        //save the selected ones
        std::vector<bool> SelV(mesh.vert.size(),false);
        for (size_t i=0;i<mesh.vert.size();i++)
            if (mesh.vert[i].IsS())SelV[i]=true;

        //update
        vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
        vcg::tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(mesh);

        //then make sure that edge sel consistent alogn adjacent faces
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!mesh.face[i].IsFaceEdgeS(j))continue;
                if (vcg::face::IsBorder(mesh.face[i],j))continue;
                FaceType *fOpp=mesh.face[i].FFp(j);
                int IOpp=mesh.face[i].FFi(j);
                assert(fOpp!=&mesh.face[i]);
                assert(fOpp!=NULL);
                assert(IOpp>=0);
                assert(IOpp<3);
                assert(fOpp->IsFaceEdgeS(IOpp));
            }

        //create star meshes
        std::vector<MeshType*> StarMeshes;
        std::vector<std::vector<int> > OrigIndex;
        std::vector<std::vector<PosType> > ProjBasis;
        GetStarMeshes(mesh,StarMeshes,ProjBasis,OrigIndex);

        //and finally smooth
        for (size_t i=0;i<UVP.Steps;i++)
            SmoothStep(mesh,StarMeshes,OrigIndex,ProjBasis,UVP,SelV);


        for (size_t i=0;i<StarMeshes.size();i++)
        {
            if (StarMeshes[i]==NULL)continue;
            delete(StarMeshes[i]);
        }

        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (SelV[i])
                mesh.vert[i].SetS();
        }
    }
}; //end

template <class MeshType>
void SmoothPolygonalMeshByLocalUV(MeshType &mesh,
                                  const std::vector<std::pair<size_t,size_t> > &SharpFeatures,
                                  const std::vector<size_t> &Corners,
                                  const size_t Steps=100)
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;

    //save positions
    std::set<std::pair<CoordType,CoordType> > SharpFPos;
    std::set<CoordType> SharpPos;
    for (size_t i=0;i<SharpFeatures.size();i++)
    {
        size_t IndexF=SharpFeatures[i].first;
        size_t IndexE=SharpFeatures[i].second;
        size_t numV=mesh.face[IndexF].VN();
        CoordType P0=mesh.face[IndexF].V(IndexE)->P();
        CoordType P1=mesh.face[IndexF].V((IndexE+1)%numV)->P();
        std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
        SharpFPos.insert(Key);
    }

    for (size_t i=0;i<Corners.size();i++)
        SharpPos.insert(mesh.vert[Corners[i]].P());

    //triangulate
    MeshType swapM;
    vcg::tri::Append<MeshType,MeshType>::Mesh(swapM,mesh);

    vcg::PolygonalAlgorithm<MeshType>::Triangulate(swapM);
    vcg::tri::UpdateTopology<MeshType>::FaceFace(swapM);
    vcg::tri::UpdateFlags<MeshType>::FaceBorderFromFF(swapM);
    vcg::tri::UpdateFlags<MeshType>::VertexBorderFromFaceBorder(swapM);

    vcg::tri::UpdateFlags<MeshType>::VertexClearS(swapM);
    vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(swapM);
    //then set it back on the mesh
    for (size_t i=0;i<swapM.face.size();i++)
    {
        size_t numV=swapM.face[i].VN();
        assert(numV==3);
        for (size_t j=0;j<3;j++)
        {
            CoordType P0=swapM.face[i].V0(j)->P();
            CoordType P1=swapM.face[i].V1(j)->P();
            std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
            if (SharpFPos.count(Key)==0)continue;
            swapM.face[i].SetFaceEdgeS(j);
        }
    }

    for (size_t i=0;i<swapM.vert.size();i++)
    {
        CoordType P=swapM.vert[i].P();
        if (SharpPos.count(P)==0)continue;
        swapM.vert[i].SetS();
    }

    //copy back values
    typename Local_Param_Smooth<MeshType>::UVSmoothParam UVP;
    UVP.LineEdgeSel=true;
    UVP.FixSel=true;
    UVP.Steps=Steps;
    Local_Param_Smooth<MeshType>::Smooth(swapM,UVP);

    //then copy bback values
    for (size_t i=0;i<mesh.vert.size();i++)
        mesh.vert[i].P()=swapM.vert[i].P();

}

#endif // LOCAL_PARAM_SMOOTH
