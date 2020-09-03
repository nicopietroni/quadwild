/****************************************************************************
**
** Copyright (C) 2011 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
**     the names of its contributors may be used to endorse or promote
**     products derived from this software without specific prior written
**     permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
** $QT_END_LICENSE$
**
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
#include <wrap/gl/trimesh.h>

#include <wrap/io_trimesh/export_ply.h>
//#include <vcg/complex/algorithms/polygonal_algorithms.h>
//#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
//#include <vcg/complex/algorithms/curve_on_manifold.h>

using namespace vcg;
class CFace;
class CVertex;

struct MyUsedTypes : public UsedTypes<	Use<CVertex>::AsVertexType,
        Use<CFace>::AsFaceType>{};

//compositing wanted proprieties
class CVertex : public vcg::Vertex< MyUsedTypes,
        vcg::vertex::TexCoord2d,
        vcg::vertex::Coord3d,
        vcg::vertex::Normal3d,
        vcg::vertex::BitFlags,
        vcg::vertex::VFAdj,
        vcg::vertex::Qualityd,
        vcg::vertex::Color4b,
        vcg::vertex::Mark,
        vcg::vertex::CurvatureDird>{};

class CFace   : public vcg::Face<  MyUsedTypes,
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

};

class CMesh   : public vcg::tri::TriMesh< std::vector<CVertex>,std::vector<CFace> >
{
public:
    std::vector<std::pair<size_t,size_t> > SharpFeatures;
    std::vector<size_t> SharpCorners;
    //ScalarType FlatDegree;

    void UpdateAttributes()
    {
        vcg::tri::UpdateNormal<CMesh>::PerFaceNormalized(*this);
        vcg::tri::UpdateNormal<CMesh>::PerVertexNormalized(*this);
        vcg::tri::UpdateBounding<CMesh>::Box(*this);
        vcg::tri::UpdateTopology<CMesh>::FaceFace(*this);
        vcg::tri::UpdateTopology<CMesh>::VertexFace(*this);
        vcg::tri::UpdateFlags<CMesh>::FaceBorderFromFF(*this);
        vcg::tri::UpdateFlags<CMesh>::VertexBorderFromFaceBorder(*this);
    }

    bool LoadField(std::string field_filename)
    {
        int position0=field_filename.find(".ffield");
        int position1=field_filename.find(".rosy");


        if (position0!=-1)
        {
            bool loaded=vcg::tri::io::ImporterFIELD<CMesh>::LoadFFIELD(*this,field_filename.c_str());
            if (!loaded)return false;
            vcg::tri::CrossField<CMesh>::OrientDirectionFaceCoherently(*this);
            vcg::tri::CrossField<CMesh>::UpdateSingularByCross(*this);
            return true;
        }
        if (position1!=-1)
        {
            std::cout<<"Importing ROSY field"<<std::endl;
            bool loaded=vcg::tri::io::ImporterFIELD<CMesh>::Load4ROSY(*this,field_filename.c_str());
            std::cout<<"Imported ROSY field"<<std::endl;
            if (!loaded)return false;
            vcg::tri::CrossField<CMesh>::OrientDirectionFaceCoherently(*this);
            vcg::tri::CrossField<CMesh>::UpdateSingularByCross(*this);
            return true;
        }
        return false;
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
            int err=vcg::tri::io::ImporterPLY<CMesh>::Open(*this,filename.c_str());
            if (err!=vcg::ply::E_NOERROR)return false;
            return true;
        }
        if (position1!=-1)
        {
            int mask;
            vcg::tri::io::ImporterOBJ<CMesh>::LoadMask(filename.c_str(),mask);
            int err=vcg::tri::io::ImporterOBJ<CMesh>::Open(*this,filename.c_str(),mask);
            if ((err!=0)&&(err!=5))return false;
            return true;
        }
        if (position2!=-1)
        {
            int err=vcg::tri::io::ImporterOFF<CMesh>::Open(*this,filename.c_str());
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
        std::cout<<"There are "<<SharpFeatures.size()<<" Sharp edges"<<std::endl;
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

    void SelectSharpFeatures()
    {
        vcg::tri::UpdateFlags<CMesh>::FaceClearFaceEdgeS(*this);
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
        fprintf(f,"%d\n",SharpCorners.size());
        for (size_t i=0;i<SharpCorners.size();i++)
            fprintf(f,"%d\n",SharpCorners[i]);
        fclose(f);
        return true;
    }

    bool SaveFeatures(std::string &FeaturePath)
    {
        FILE *f=NULL;
        f=fopen(FeaturePath.c_str(),"wt");
        if(f==NULL) return false;
        fprintf(f,"%d\n",SharpFeatures.size());
        for (size_t i=0;i<SharpFeatures.size();i++)
            fprintf(f,"%d,%d\n",SharpFeatures[i].first,SharpFeatures[i].second);
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

//    int GenusOfSelectedFaces()
//    {
//        vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(*this);
//        std::set<std::pair<size_t,size_t> > EdgeSet;
//        size_t NumF=0;
//        size_t NumV=0;
//        size_t NumE=0;
//        for (size_t i=0;i<face.size();i++)
//        {
//            if (!face[i].IsS())continue;
//            NumF++;
//            for (size_t j=0;j<3;j++)
//            {
//               size_t IndV0=vcg::tri::Index(*this,face[i].V0(j));
//               size_t IndV1=vcg::tri::Index(*this,face[i].V1(j));
//               EdgeSet.insert(std::pair<size_t,size_t>(std::min(IndV0,IndV1),std::max(IndV0,IndV1)));
//            }
//        }
//        for (size_t i=0;i<vert.size();i++)
//        {
//             if (vert[i].IsD())continue;
//            if (vert[i].IsS())NumV++;
//        }

//        NumE=EdgeSet.size();
//        return ( NumV + NumF - NumE );
//    }

    void SelectPos(const  std::vector<std::vector<vcg::face::Pos<FaceType> > > &ToSel,bool SetSel)
    {
        for (size_t i=0;i<ToSel.size();i++)
           SelectPos(ToSel[i],SetSel);
    }

    void GLDrawSharpEdges()
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9999);
        glLineWidth(10);
        glBegin(GL_LINES);
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;

                vcg::glColor(vcg::Color4b(255,0,0,255));

                CoordType Pos0=face[i].P0(j);
                CoordType Pos1=face[i].P1(j);
                vcg::glVertex(Pos0);
                vcg::glVertex(Pos1);
            }
        glEnd();
        glPopAttrib();
    }
};

#endif
