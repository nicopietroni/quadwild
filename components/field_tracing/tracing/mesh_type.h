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

#ifndef MESH_TYPE_H
#define MESH_TYPE_H

/// vcg imports
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/refine.h>
/// wrapper imports
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/import_field.h>
//#include <wrap/gl/trimesh.h>

#include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
#include <vcg/complex/algorithms/parametrization/uv_utils.h>
//#include <vcg/complex/algorithms/polygonal_algorithms.h>
//#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
//#include <vcg/complex/algorithms/curve_on_manifold.h>

using namespace vcg;
class TraceFace;
class TraceVertex;

struct MyUsedTypes : public UsedTypes<	Use<TraceVertex>::AsVertexType,
        Use<TraceFace>::AsFaceType>{};

//compositing wanted proprieties
class TraceVertex : public vcg::Vertex< MyUsedTypes,
        vcg::vertex::TexCoord2d,
        vcg::vertex::Coord3d,
        vcg::vertex::Normal3d,
        vcg::vertex::BitFlags,
        vcg::vertex::VFAdj,
        vcg::vertex::Qualityd,
        vcg::vertex::Color4b,
        vcg::vertex::Mark,
        vcg::vertex::CurvatureDird>
{
public:
#ifdef MULTI_FRAME
    std::vector<CoordType> FramePos;
#endif
    CoordType RPos;
    size_t SingularityValence;

    void ImportData(const TraceVertex  & left )
    {
        vcg::Vertex< MyUsedTypes,
                vcg::vertex::TexCoord2d,
                vcg::vertex::Coord3d,
                vcg::vertex::Normal3d,
                vcg::vertex::BitFlags,
                vcg::vertex::VFAdj,
                vcg::vertex::Qualityd,
                vcg::vertex::Color4b,
                vcg::vertex::Mark,
                vcg::vertex::CurvatureDird>::ImportData(left);

        RPos=left.RPos;
        SingularityValence=left.SingularityValence;

#ifdef MULTI_FRAME
    FramePos=left.FramePos;
#endif
    }

};

class TraceFace   : public vcg::Face<  MyUsedTypes,
        vcg::face::VertexRef,
        vcg::face::Normal3d,
        vcg::face::BitFlags,
        vcg::face::CurvatureDird,
        vcg::face::FFAdj,
        vcg::face::VFAdj,
        vcg::face::Qualityd,
        vcg::face::Color4b,
        vcg::face::Mark,
        vcg::face::CurvatureDird,
        vcg::face::WedgeTexCoord2d>
{
public:
    bool FullTraced;
    size_t IndexOriginal;

    void ImportData(const TraceFace  & left )
    {
        vcg::Face<  MyUsedTypes,
                vcg::face::VertexRef,
                vcg::face::Normal3d,
                vcg::face::BitFlags,
                vcg::face::CurvatureDird,
                vcg::face::FFAdj,
                vcg::face::VFAdj,
                vcg::face::Qualityd,
                vcg::face::Color4b,
                vcg::face::Mark,
                vcg::face::CurvatureDird,
                vcg::face::WedgeTexCoord2d>::ImportData(left);

        FullTraced=left.FullTraced;
        IndexOriginal=left.IndexOriginal;
    }
};

class TraceMesh   : public vcg::tri::TriMesh< std::vector<TraceVertex>,std::vector<TraceFace> >
{
public:
    std::vector<std::pair<size_t,size_t> > SharpFeatures;
    std::vector<size_t> SharpCorners;


    void UpdateAttributes()
    {
        vcg::tri::UpdateNormal<TraceMesh>::PerFaceNormalized(*this);
        vcg::tri::UpdateNormal<TraceMesh>::PerVertexNormalized(*this);
        vcg::tri::UpdateBounding<TraceMesh>::Box(*this);
        vcg::tri::UpdateTopology<TraceMesh>::FaceFace(*this);
        vcg::tri::UpdateTopology<TraceMesh>::VertexFace(*this);
        vcg::tri::UpdateFlags<TraceMesh>::FaceBorderFromFF(*this);
        vcg::tri::UpdateFlags<TraceMesh>::VertexBorderFromFaceBorder(*this);
    }

    bool LoadField(std::string field_filename)
    {
        int position0=field_filename.find(".ffield");
        int position1=field_filename.find(".rosy");


        if (position0!=-1)
        {
            bool loaded=vcg::tri::io::ImporterFIELD<TraceMesh>::LoadFFIELD(*this,field_filename.c_str());
            if (!loaded)return false;
            vcg::tri::CrossField<TraceMesh>::OrientDirectionFaceCoherently(*this);
            vcg::tri::CrossField<TraceMesh>::UpdateSingularByCross(*this,true);
            return true;
        }
        if (position1!=-1)
        {
            std::cout<<"Importing ROSY field"<<std::endl;
            bool loaded=vcg::tri::io::ImporterFIELD<TraceMesh>::Load4ROSY(*this,field_filename.c_str());
            std::cout<<"Imported ROSY field"<<std::endl;
            if (!loaded)return false;
            vcg::tri::CrossField<TraceMesh>::OrientDirectionFaceCoherently(*this);
            vcg::tri::CrossField<TraceMesh>::UpdateSingularByCross(*this,true);
            return true;
        }
        return false;
    }

    void ScatterColorByQualityFace()
    {
        int MaxV=0;
        for (size_t i=0;i<face.size();i++)
            MaxV=std::max(MaxV,(int)face[i].Q());

        for (size_t i=0;i<face.size();i++)
        {
            vcg::Color4b CurrCol=vcg::Color4b::Scatter((int)MaxV+1,(int)face[i].Q());
            face[i].C()=CurrCol;
        }
    }

    bool LoadMesh(std::string filename)
    {
        Clear();
        if(filename.empty()) return false;
        int position0=filename.find(".ply");
        int position1=filename.find(".obj");
        int position2=filename.find(".off");

        if (position0!=-1)
        {
            int err=vcg::tri::io::ImporterPLY<TraceMesh>::Open(*this,filename.c_str());
            if (err!=vcg::ply::E_NOERROR)return false;
            return true;
        }
        if (position1!=-1)
        {
            int mask;
            vcg::tri::io::ImporterOBJ<TraceMesh>::LoadMask(filename.c_str(),mask);
            int err=vcg::tri::io::ImporterOBJ<TraceMesh>::Open(*this,filename.c_str(),mask);
            if ((err!=0)&&(err!=5))return false;
            return true;
        }
        if (position2!=-1)
        {
            int err=vcg::tri::io::ImporterOFF<TraceMesh>::Open(*this,filename.c_str());
            if (err!=0)return false;
            return true;
        }
        return false;
    }
    bool SaveTriMesh(const std::string &filename)
    {
        if(filename.empty()) return false;
        int position0=filename.find(".ply");
        int position1=filename.find(".obj");
        int position2=filename.find(".off");

        if (position0!=-1)
        {
            int mask=vcg::tri::io::Mask::IOM_FACECOLOR|vcg::tri::io::Mask::IOM_WEDGTEXCOORD;
            int err=vcg::tri::io::ExporterPLY<TraceMesh>::Save(*this,filename.c_str(),mask);
            if (err!=vcg::ply::E_NOERROR)return false;
            return true;
        }
        if (position1!=-1)
        {
            int mask=vcg::tri::io::Mask::IOM_FACECOLOR|vcg::tri::io::Mask::IOM_WEDGTEXCOORD;
            int err=vcg::tri::io::ExporterOBJ<TraceMesh>::Save(*this,filename.c_str(),mask);
            if ((err!=0)&&(err!=5))return false;
            return true;
        }
        if (position2!=-1)
        {
            int err=vcg::tri::io::ExporterOFF<TraceMesh>::Save(*this,filename.c_str());
            if (err!=0)return false;
            return true;
        }
        return false;
    }

    void UpdateFromCoordPairs(const std::vector<std::pair<CoordType,CoordType> > &SharpCoords,
                              bool check=true)
    {
        std::map<std::pair<CoordType,CoordType>,std::pair<size_t,size_t> > EdgeMap;
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                CoordType P0=face[i].P0(j);
                CoordType P1=face[i].P1(j);
                std::pair<CoordType,CoordType> key(std::min(P0,P1),std::max(P0,P1));
                EdgeMap[key]=std::pair<size_t,size_t>(i,j);
            }
        SharpFeatures.clear();
        for (size_t i=0;i<SharpCoords.size();i++)
        {
            CoordType P0=SharpCoords[i].first;
            CoordType P1=SharpCoords[i].second;
            std::pair<CoordType,CoordType> key(std::min(P0,P1),std::max(P0,P1));
            if (check)
            {
                assert(EdgeMap.count(key)>0);
                SharpFeatures.push_back(EdgeMap[key]);
            }
            else
            {
                if (EdgeMap.count(key)>0)
                    SharpFeatures.push_back(EdgeMap[key]);
            }
        }
        //std::cout<<"There are "<<SharpFeatures.size()<<" Sharp edges"<<std::endl;
    }

    void GetSharpCoordPairs(std::vector<std::pair<CoordType,CoordType> > &SharpCoords)
    {
        SharpCoords.clear();
        for (size_t i=0;i<SharpFeatures.size();i++)
        {
            size_t IndexF=SharpFeatures[i].first;
            size_t IndexE=SharpFeatures[i].second;
            CoordType P0=face[IndexF].P0(IndexE);
            CoordType P1=face[IndexF].P1(IndexE);
            std::pair<CoordType,CoordType> key(std::min(P0,P1),std::max(P0,P1));
            SharpCoords.push_back(key);
        }
        std::sort(SharpCoords.begin(),SharpCoords.end());
        auto last=std::unique(SharpCoords.begin(),SharpCoords.end());
        SharpCoords.erase(last, SharpCoords.end());
    }

    size_t NumDuplicatedV()
    {
        std::set<CoordType> Pos;
        vcg::tri::UpdateSelection<MeshType>::VertexClear(*this);
        size_t numDupl=0;
        for (size_t i=0;i<vert.size();i++)
        {
            if (Pos.count(vert[i].P())>0)
            {
                vert[i].SetS();
                numDupl++;
            }
            Pos.insert(vert[i].P());
        }
        return numDupl;
    }

    void Perturb(VertexType &v,ScalarType Magnitudo)
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
        v.RPos=v.P();
    }

    bool RepositionDuplicatedV()
    {
        size_t NumD=NumDuplicatedV();
        if (NumD==0)return false;
        //int dilate_step=0;
        ScalarType Magnitudo=2;
        do
        {
            std::cout<<"Repositioning "<<NumD<<" duplicated vertices"<<std::endl;


            //            dilate_step++;
            //            for (size_t i=0;i<dilate_step;i++)
            //            {
            //                vcg::tri::UpdateSelection<MeshType>::FaceFromVertexLoose(*this);
            //                vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(*this);
            //            }
            for (size_t i=0;i<vert.size();i++)
                if (vert[i].IsS())Perturb(vert[i],Magnitudo);

            Magnitudo*=2;
            //vcg::tri::Smooth<MeshType>::VertexCoordLaplacian(*this,1,true);
            NumD=NumDuplicatedV();
        }
        while(NumD>0);
        vcg::tri::UpdateBounding<MeshType>::Box(*this);
        vcg::tri::UpdateNormal<MeshType>::PerVertexNormalizedPerFace(*this);
        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(*this);
        return true;
    }

    bool RemoveZeroAreaF()
    {
        //        int nonManifV=0;
        //        int degF=0;

        int zeroAFace=0;
        bool modified=false;
        ScalarType Magnitudo=2;
        do{
            modified=false;
            for (size_t i=0;i<face.size();i++)
            {
                if (vcg::DoubleArea(face[i])>0)continue;
                Perturb(*face[i].V(0),Magnitudo);
                Perturb(*face[i].V(1),Magnitudo);
                Perturb(*face[i].V(2),Magnitudo);
                modified=true;
                zeroAFace++;
            }
            Magnitudo*=2;
        }while (modified);
        vcg::tri::Allocator<MeshType>::CompactEveryVector(*this);
        std::cout<<"Adjusted "<<zeroAFace<<" zero area faces"<<std::endl;
        //        std::cout<<"Removed "<<degF<<" degenerate faces"<<std::endl;
        //        std::cout<<"Removed "<<zeroAFace<<" nonManifV "<<std::endl;
        UpdateAttributes();
        return modified;
    }


    //    bool RemoveZeroAreaF()
    //    {
    //        int nonManifV=0;
    //        int degF=0;
    //        int zeroAFace=0;
    //        bool modified=false;
    //        do{
    //            zeroAFace=vcg::tri::Clean<MeshType>::RemoveZeroAreaFace(*this);
    //            degF=vcg::tri::Clean<MeshType>::RemoveDegenerateFace(*this);
    //            nonManifV=vcg::tri::Clean<MeshType>::RemoveNonManifoldVertex(*this);
    //            vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(*this);
    //            modified|=((zeroAFace!=0)&&(degF!=0)&&(nonManifV!=0));
    //        }while ((zeroAFace!=0)&&(degF!=0)&&(nonManifV!=0));
    //        vcg::tri::Allocator<MeshType>::CompactEveryVector(*this);
    //        std::cout<<"Removed "<<zeroAFace<<" zero area faces"<<std::endl;
    //        std::cout<<"Removed "<<degF<<" degenerate faces"<<std::endl;
    //        std::cout<<"Removed "<<zeroAFace<<" nonManifV "<<std::endl;
    //        UpdateAttributes();
    //        return modified;
    //    }

    void SolveGeometricIssues()
    {
        srand(0);
        bool HasRepositioned=false;
        bool HasSolvedZeroF=false;
        do{
            HasRepositioned=RepositionDuplicatedV();
            HasSolvedZeroF=RemoveZeroAreaF();
        }while (HasRepositioned || HasSolvedZeroF);
        UpdateAttributes();

    }

    void UpdateSharpFeaturesFromSelection()
    {
        SharpFeatures.clear();
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                SharpFeatures.push_back(std::pair<size_t,size_t>(i,j));
                //std::cout<<"DE BOIA"<<std::endl;
            }
    }

    void SelectSharpFeatures()
    {
        vcg::tri::UpdateFlags<TraceMesh>::FaceClearFaceEdgeS(*this);
        for (size_t i=0;i<SharpFeatures.size();i++)
        {
            int IndexF=SharpFeatures[i].first;
            int IndexE=SharpFeatures[i].second;
            face[IndexF].SetFaceEdgeS(IndexE);
            if (vcg::face::IsBorder(face[IndexF],IndexE))continue;
            FaceType *Fopp=face[IndexF].FFp(IndexE);
            int IndexOpp=face[IndexF].FFi(IndexE);
            Fopp->SetFaceEdgeS(IndexOpp);
        }
    }

    bool LoadSharpFeatures(std::string &FeaturePath)
    {
        SharpFeatures.clear();
        FILE *f=NULL;
        f=fopen(FeaturePath.c_str(),"rt");
        if(f==NULL) return false;
        int Num=0;
        fscanf(f,"%d\n",&Num);
        std::cout<<"Num "<<Num<<std::endl;
        for (size_t i=0;i<(size_t)Num;i++)
        {
            int FIndex,EIndex;
            int FType;
            fscanf(f,"%d,%d,%d\n",&FType,&FIndex,&EIndex);
            assert(FIndex>=0);
            assert(FIndex<(int)face.size());
            assert(EIndex>=0);
            assert(EIndex<4);
            face[FIndex].SetFaceEdgeS(EIndex);
            SharpFeatures.push_back(std::pair<size_t,size_t>(FIndex,EIndex));
        }
        fclose(f);
        return true;
    }

    bool SaveSharpCorners(std::string &SharpCornerPath)
    {
        FILE *f=NULL;
        f=fopen(SharpCornerPath.c_str(),"wt");
        if(f==NULL) return false;
        fprintf(f,"%d\n",(int)SharpCorners.size());
        for (size_t i=0;i<SharpCorners.size();i++)
            fprintf(f,"%d\n",(int)SharpCorners[i]);
        fclose(f);
        return true;
    }

    bool SaveFeatures(std::string &FeaturePath)
    {
        FILE *f=NULL;
        f=fopen(FeaturePath.c_str(),"wt");
        if(f==NULL) return false;
        fprintf(f,"%d\n",(int)SharpFeatures.size());
        for (size_t i=0;i<SharpFeatures.size();i++)
            fprintf(f,"%d,%d\n",(int)SharpFeatures[i].first,(int)SharpFeatures[i].second);
        fclose(f);
        return true;
    }

    void SelectPos(const  std::vector<vcg::face::Pos<FaceType> > &ToSel,bool SetSel)
    {
        for (size_t i=0;i<ToSel.size();i++)
        {
            vcg::face::Pos<FaceType> Pos0=ToSel[i];
            if (SetSel)
            {
                Pos0.F()->SetFaceEdgeS(Pos0.E());
                Pos0.FlipF();
                Pos0.F()->SetFaceEdgeS(Pos0.E());
            }
            else
            {
                Pos0.F()->ClearFaceEdgeS(Pos0.E());
                Pos0.FlipF();
                Pos0.F()->ClearFaceEdgeS(Pos0.E());
            }
        }
    }


    void SelectPos(const  std::vector<std::vector<vcg::face::Pos<FaceType> > > &ToSel,bool SetSel)
    {
        for (size_t i=0;i<ToSel.size();i++)
            SelectPos(ToSel[i],SetSel);
    }

    //    void GLDrawSharpEdges(vcg::Color4b col=vcg::Color4b(255,0,0,255),
    //                          ScalarType GLSize=5)
    //    {
    //        glPushAttrib(GL_ALL_ATTRIB_BITS);
    //        glDisable(GL_LIGHTING);
    //        glDepthRange(0,0.9999);
    //        glLineWidth(GLSize);
    //        vcg::glColor(col);
    //        glBegin(GL_LINES);
    //        for (size_t i=0;i<SharpFeatures.size();i++)
    //        {

    //            size_t IndexF=SharpFeatures[i].first;
    //            size_t IndexE=SharpFeatures[i].second;
    //            CoordType Pos0=face[IndexF].P0(IndexE);
    //            CoordType Pos1=face[IndexF].P1(IndexE);
    //            vcg::glVertex(Pos0);
    //            vcg::glVertex(Pos1);
    //        }
    //        glEnd();
    //        glPopAttrib();
    //    }

    void InitSingVert()
    {
        // query if an attribute is present or not
        bool hasSingular = vcg::tri::HasPerVertexAttribute(*this,std::string("Singular"));
        bool hasSingularIndex = vcg::tri::HasPerVertexAttribute(*this,std::string("SingularIndex"));

        assert(hasSingular);
        assert(hasSingularIndex);

        typename MeshType::template PerVertexAttributeHandle<bool> Handle_Singular;
        Handle_Singular=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<bool>(*this,std::string("Singular"));
        typename MeshType::template PerVertexAttributeHandle<int> Handle_SingularIndex;
        Handle_SingularIndex =vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(*this,std::string("SingularIndex"));

        for (size_t i=0;i<vert.size();i++)
        {
            vert[i].SingularityValence=4;
            if (vert[i].IsD())continue;
            if (!Handle_Singular[i])continue;

            int SingIndex=Handle_SingularIndex[i];

            switch (SingIndex)
            {
            case 1:vert[i].SingularityValence=5;break;
            case 2:vert[i].SingularityValence=6;break;
            case 3:vert[i].SingularityValence=3;break;
            case 4:vert[i].SingularityValence=2;break;
            default:break;
            }
        }
    }

    void WichFaceEdge(const size_t &IndexV0,
                      const size_t &IndexV1,
                      int &IndexF,
                      int &IndexE)
    {
        IndexF=-1;
        IndexE=-1;
        std::pair<size_t,size_t> targetE=std::pair<size_t,size_t>(std::min(IndexV0,IndexV1),
                                                                  std::max(IndexV0,IndexV1));
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                size_t TestV0=vcg::tri::Index(*this,face[i].V0(j));
                size_t TestV1=vcg::tri::Index(*this,face[i].V1(j));
                std::pair<size_t,size_t> testE=std::pair<size_t,size_t>(std::min(TestV0,TestV1),
                                                                        std::max(TestV0,TestV1));
                if (targetE==testE)
                {
                    IndexF=i;
                    IndexE=j;
                    return;
                }
            }
    }

    CoordType MoveCenterOnZero()
    {
        CoordType Center=bbox.Center();
        for (size_t i=0;i<vert.size();i++)
            vert[i].P()-=Center;
        UpdateAttributes();
        return Center;
    }

    void RestoreRPos()
    {
        for (size_t i=0;i<vert.size();i++)
            vert[i].P()=vert[i].RPos;
    }

    void InitRPos()
    {
        for (size_t i=0;i<vert.size();i++)
            vert[i].RPos=vert[i].P();
    }


    //    CoordType UVTo3DPos(const vcg::Point2<ScalarType> &UVPos)
    //    {
    //        CoordType Pos(0,0,0);
    //        Pos.X()=UVPos.X();
    //        Pos.Y()=UVPos.Y();
    //        return Pos;
    //    }

    //    void GLDrawEdgeUV()
    //    {
    //        glPushAttrib(GL_ALL_ATTRIB_BITS);
    //        glPushMatrix();
    //        vcg::Box2<ScalarType> uv_box=vcg::tri::UV_Utils<TraceMesh>::PerWedgeUVBox(*this);
    //        vcg::glScale(3.0f/uv_box.Diag());
    //        vcg::glTranslate(CoordType(-uv_box.Center().X(),
    //                                   -uv_box.Center().Y(),0));
    //        glDisable(GL_LIGHTING);
    //        glDisable(GL_LIGHT0);
    //        glLineWidth(10);
    //        //glDepthRange(0,0.9999);
    //        glBegin(GL_LINES);
    //        for (size_t i=0;i<face.size();i++)
    //        {
    //            vcg::glColor(face[i].C());
    //            CoordType Pos[3];
    //            Pos[0]=UVTo3DPos(face[i].WT(0).P());
    //            Pos[1]=UVTo3DPos(face[i].WT(1).P());
    //            Pos[2]=UVTo3DPos(face[i].WT(2).P());
    //            for (size_t j=0;j<3;j++)
    //            {
    //                if (!face[i].IsFaceEdgeS(j))continue;
    //                vcg::glColor(vcg::Color4b(0,0,0,255));
    //                vcg::glVertex(Pos[j]);
    //                vcg::glVertex(Pos[(j+1)%3]);
    //            }
    //        }
    //        glEnd();
    //        glPopMatrix();
    //        glPopAttrib();
    //    }

    //    void GLDrawUV(int TxtIndex=-1,bool colorPerVert=false)
    //    {
    ////        vcg::Box2<ScalarType> uv_box=vcg::tri::UV_Utils<TraceMesh>::PerWedgeUVBox(deformed_mesh);
    ////        ScalarType UVScale=3.0f/uv_box.Diag();
    ////        for (size_t i=0;i<deformed_mesh.face.size();i++)
    ////        {
    ////            deformed_mesh.face[i].WT(0).P()*=UVScale;
    ////            deformed_mesh.face[i].WT(1).P()*=UVScale;
    ////            deformed_mesh.face[i].WT(2).P()*=UVScale;
    ////        }

    //        glPushAttrib(GL_ALL_ATTRIB_BITS);
    //        glPushMatrix();
    //        vcg::Box2<ScalarType> uv_box=vcg::tri::UV_Utils<TraceMesh>::PerWedgeUVBox(*this);
    //        vcg::glScale(3.0f/uv_box.Diag());
    //        //ScalarType UVScale=3.0f/uv_box.Diag();
    //        vcg::glTranslate(CoordType(-uv_box.Center().X(),
    //                                   -uv_box.Center().Y(),0));
    //        glDisable(GL_LIGHTING);
    //        glDisable(GL_LIGHT0);

    //        if (TxtIndex>=0)
    //        {
    //            glActiveTexture(GL_TEXTURE0);
    //            glEnable(GL_TEXTURE_2D);
    //            glBindTexture(GL_TEXTURE_2D, TxtIndex);
    //        }
    //        glBegin(GL_TRIANGLES);
    //        for (size_t i=0;i<face.size();i++)
    //        {
    //            if (!colorPerVert)
    //                vcg::glColor(face[i].C());
    //            CoordType Pos0=UVTo3DPos(face[i].WT(0).P());
    //            CoordType Pos1=UVTo3DPos(face[i].WT(1).P());
    //            CoordType Pos2=UVTo3DPos(face[i].WT(2).P());

    //            if (TxtIndex>=0)
    //                vcg::glTexCoord(face[i].WT(0).P());
    //            if (colorPerVert)
    //                vcg::glColor(face[i].V(0)->C());

    //            vcg::glVertex(Pos0);
    //            if (TxtIndex>=0)
    //                vcg::glTexCoord(face[i].WT(1).P());
    //            if (colorPerVert)
    //                vcg::glColor(face[i].V(1)->C());

    //            vcg::glVertex(Pos1);
    //            if (TxtIndex>=0)
    //                vcg::glTexCoord(face[i].WT(2).P());
    //            if (colorPerVert)
    //                vcg::glColor(face[i].V(2)->C());

    //            vcg::glVertex(Pos2);
    //        }
    //        glEnd();
    //        glPopMatrix();
    //        glPopAttrib();
    //    }


    bool SaveOrigFace(const std::string &filename)
    {
        if(filename.empty()) return false;
        std::ofstream myfile;
        myfile.open (filename.c_str());

        for (size_t i=0;i<face.size();i++)
        {
            myfile <<face[i].IndexOriginal <<std::endl;
        }

        myfile.close();
        return true;
    }


    bool LoadOrigFaces(std::string &filename)
    {
        FILE *f=NULL;
        f=fopen(filename.c_str(),"rt");
        if(f==NULL) return false;
        for (size_t i=0;i<face.size();i++)
        {
            int FIndex;
            fscanf(f,"%d\n",&FIndex);
            assert(i>=0);
            assert(i<face.size());
            face[i].IndexOriginal=(size_t)FIndex;
        }
        fclose(f);
        return true;
    }
};

#endif
