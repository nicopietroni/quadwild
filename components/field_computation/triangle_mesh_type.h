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

#ifndef MY_TRI_MESH_TYPE
#define MY_TRI_MESH_TYPE

//#define MINDOT -0.99

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/simplex/face/topology.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
//#include <wrap/gl/trimesh.h>
//
#include <wrap/io_trimesh/export_field.h>
#include <wrap/io_trimesh/import_field.h>
#include <iostream>
#include <fstream>
#include <vcg/complex/algorithms/attribute_seam.h>
#include <vcg/complex/algorithms/crease_cut.h>
#include "fields/field_smoother.h"


class FieldTriFace;
//class MyTriEdge;
class FieldTriVertex;

enum FeatureKind{ETConcave,ETConvex,ETNone};

struct TriUsedTypes: public vcg::UsedTypes<vcg::Use<FieldTriVertex>::AsVertexType,
        //vcg::Use<MyTriEdge>::AsEdgeType,
        vcg::Use<FieldTriFace>::AsFaceType>{};

class FieldTriVertex:public vcg::Vertex<TriUsedTypes,
        vcg::vertex::Coord3d,
        vcg::vertex::Color4b,
        vcg::vertex::Normal3d,
        vcg::vertex::VFAdj,
        vcg::vertex::BitFlags,
        vcg::vertex::CurvatureDird,
        vcg::vertex::Qualityd,
        vcg::vertex::TexCoord2d,
        vcg::vertex::Mark>
{
};

class FieldTriFace:public vcg::Face<TriUsedTypes,
        vcg::face::VertexRef,
        vcg::face::VFAdj,
        vcg::face::FFAdj,
        vcg::face::BitFlags,
        vcg::face::Normal3d,
        vcg::face::CurvatureDird,
        vcg::face::Color4b,
        vcg::face::Qualityd,
        vcg::face::WedgeTexCoord2d,
        vcg::face::Mark>
{

public:

    size_t IndexOriginal;

    void ImportData(const FieldTriFace  & left )
    {
        IndexOriginal=left.IndexOriginal;
        vcg::Face<TriUsedTypes,
                vcg::face::VertexRef,
                vcg::face::VFAdj,
                vcg::face::FFAdj,
                vcg::face::BitFlags,
                vcg::face::Normal3d,
                vcg::face::CurvatureDird,
                vcg::face::Color4b,
                vcg::face::Qualityd,
                vcg::face::WedgeTexCoord2d,
                vcg::face::Mark>::ImportData(left);
    }

    FeatureKind FKind[3];
};

//enum GoemPrecondition{NOVertManifold,NOFaceManifold,
//                      DegenerateFace,DegenerateVertex,
//                      UnreferencedVert,
//                      IsOk};


class FieldTriMesh: public vcg::tri::TriMesh< std::vector<FieldTriVertex>,
        //std::vector<MyTriEdge>,
        std::vector<FieldTriFace > >
{
    typedef std::pair<CoordType,CoordType> CoordPair;
    std::set< CoordPair > FeaturesCoord;

public:

    ScalarType LimitConcave;

    void InitEdgeType()
    {
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (IsConcaveEdge(face[i],j))
                    face[i].FKind[j]=ETConcave;
                else
                    face[i].FKind[j]=ETConvex;
            }
    }

    void InitFeatureCoordsTable()
    {
        FeaturesCoord.clear();
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                CoordPair PosEdge(std::min(face[i].P0(j),face[i].P1(j)),
                                  std::max(face[i].P0(j),face[i].P1(j)));
                FeaturesCoord.insert(PosEdge);
            }
        }
    }


    void SetFeatureFromTable()
    {
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                face[i].ClearFaceEdgeS(j);
                CoordPair PosEdge(std::min(face[i].P0(j),face[i].P1(j)),
                                  std::max(face[i].P0(j),face[i].P1(j)));
                if(FeaturesCoord.count(PosEdge)==0)continue;
                face[i].SetFaceEdgeS(j);
            }
        }
    }

    bool IsConcaveEdge(const FaceType &f0,int IndexE)
    {
        FaceType *f1=f0.cFFp(IndexE);
        if (f1==&f0)return false;
        CoordType N0=f0.cN();
        CoordType N1=f1->cN();
        CoordType EdgeDir=f0.cP1(IndexE)-f0.cP0(IndexE);
        EdgeDir.Normalize();
        CoordType Cross=N0^N1;
        return ((Cross*EdgeDir)<LimitConcave);
    }



public:

    bool LoadTriMesh(const std::string &filename,bool &allQuad)
    {
        allQuad=false;
        Clear();
        if(filename.empty()) return false;
        int position0=filename.find(".ply");
        int position1=filename.find(".obj");
        int position2=filename.find(".off");

        if (position0!=-1)
        {
            int err=vcg::tri::io::ImporterPLY<FieldTriMesh>::Open(*this,filename.c_str());
            if (err!=vcg::ply::E_NOERROR)return false;
            return true;
        }
        if (position1!=-1)
        {
#ifdef MIQ_QUADRANGULATE
            PMesh pmesh;
            int Mask;
            vcg::tri::io::ImporterOBJ<PMesh>::LoadMask(filename.c_str(), Mask);
            int err=vcg::tri::io::ImporterOBJ<PMesh>::Open(pmesh,filename.c_str(),Mask);
            if ((err!=0)&&(err!=5))return false;
            check if all quad
                    allQuad=true;
            for (size_t i=0;i<pmesh.face.size();i++)
            {
                if (pmesh.face[i].VN()==4)continue;
                allQuad=false;
            }
            if (allQuad)
            {
                for (size_t i=0;i<pmesh.face.size();i++)
                {
                    CoordType PD1[2];
                    PD1[0]=(pmesh.face[i].V(0)->P()-pmesh.face[i].V(1)->P());
                    PD1[1]=(pmesh.face[i].V(3)->P()-pmesh.face[i].V(2)->P());
                    PD1[0].Normalize();
                    PD1[1].Normalize();
                    CoordType PD2[2];
                    PD2[0]=(pmesh.face[i].V(2)->P()-pmesh.face[i].V(1)->P());
                    PD2[1]=(pmesh.face[i].V(3)->P()-pmesh.face[i].V(0)->P());
                    PD2[0].Normalize();
                    PD2[1].Normalize();
                    pmesh.face[i].PD1()=PD1[0]+PD1[1];
                    pmesh.face[i].PD2()=PD2[0]+PD2[1];
                    pmesh.face[i].PD1().Normalize();
                    pmesh.face[i].PD2().Normalize();
                }
                size_t size=pmesh.fn;
                //vcg::PolygonalAlgorithm<PMesh>::Triangulate(pmesh,false);

                pmesh.TriangulateQuadBySplit();
                for (size_t i=0;i<size;i++)
                {
                    pmesh.face[i+size].PD1()=pmesh.face[i].PD1();
                    pmesh.face[i+size].PD2()=pmesh.face[i].PD2();
                }
                //then copy the field
                Clear();
                vcg::tri::Allocator<FieldTriMesh>::AddVertices(*this,pmesh.vn);
                vcg::tri::Allocator<FieldTriMesh>::AddFaces(*this,pmesh.fn);

                for (size_t i=0;i<pmesh.vert.size();i++)
                    vert[i].P()=pmesh.vert[i].P();

                for (size_t i=0;i<pmesh.face.size();i++)
                {
                    size_t IndexV0=vcg::tri::Index(pmesh,pmesh.face[i].V(0));
                    size_t IndexV1=vcg::tri::Index(pmesh,pmesh.face[i].V(1));
                    size_t IndexV2=vcg::tri::Index(pmesh,pmesh.face[i].V(2));
                    face[i].V(0)=&vert[IndexV0];
                    face[i].V(1)=&vert[IndexV1];
                    face[i].V(2)=&vert[IndexV2];
                    face[i].PD1()=pmesh.face[i].PD1();
                    face[i].PD2()=pmesh.face[i].PD2();
                }
                UpdateDataStructures();
                for (size_t i=0;i<face.size();i++)
                {
                    face[i].PD1()-=face[i].N()*(face[i].PD1()*face[i].N());
                    face[i].PD2()-=face[i].N()*(face[i].PD2()*face[i].N());
                    face[i].PD1().Normalize();
                    face[i].PD2().Normalize();
                    CoordType Avg=face[i].PD1()+face[i].PD2();
                    Avg.Normalize();
                    CoordType Avg1=face[i].N()^Avg;
                    Avg1.Normalize();
                    face[i].PD1()=Avg+Avg1;
                    face[i].PD1().Normalize();
                    face[i].PD2()=face[i].N()^face[i].PD1();
                }
                vcg::tri::CrossField<FieldTriMesh>::OrientDirectionFaceCoherently(*this);
                vcg::tri::CrossField<FieldTriMesh>::UpdateSingularByCross(*this);
                return true;
            }
            else
            {
#endif
                int mask;
                vcg::tri::io::ImporterOBJ<FieldTriMesh>::LoadMask(filename.c_str(),mask);
                int err=vcg::tri::io::ImporterOBJ<FieldTriMesh>::Open(*this,filename.c_str(),mask);

                if ((err!=0)&&(err!=5))return false;
                return true;
#ifdef MIQ_QUADRANGULATE
            }
#endif

        }
        if (position2!=-1)
        {
            int err=vcg::tri::io::ImporterOFF<FieldTriMesh>::Open(*this,filename.c_str());
            if (err!=0)return false;
            return true;
        }
        return false;
    }

    bool SaveSharpFeatures(const std::string &filename)const
    {
        if(filename.empty()) return false;
        std::ofstream myfile;
        myfile.open (filename.c_str());
        size_t num=0;
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                num++;
            }
        myfile <<num<<std::endl;
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                if (face[i].FKind[j]==ETConcave)
                    myfile <<"0,"<< i <<","<<j<<std::endl;
                else
                    myfile <<"1,"<< i <<","<<j<<std::endl;
            }
        myfile.close();
        return true;
    }

    bool LoadField(std::string field_filename)
    {
        int position0=field_filename.find(".ffield");
        int position1=field_filename.find(".rosy");


        if (position0!=-1)
        {
            bool loaded=vcg::tri::io::ImporterFIELD<FieldTriMesh>::LoadFFIELD(*this,field_filename.c_str());
            if (!loaded)return false;
            vcg::tri::CrossField<FieldTriMesh>::OrientDirectionFaceCoherently(*this);
            vcg::tri::CrossField<FieldTriMesh>::UpdateSingularByCross(*this,true);
            return true;
        }
        if (position1!=-1)
        {
            std::cout<<"Importing ROSY field"<<std::endl;
            bool loaded=vcg::tri::io::ImporterFIELD<FieldTriMesh>::Load4ROSY(*this,field_filename.c_str());
            std::cout<<"Imported ROSY field"<<std::endl;
            if (!loaded)return false;
            vcg::tri::CrossField<FieldTriMesh>::OrientDirectionFaceCoherently(*this);
            vcg::tri::CrossField<FieldTriMesh>::UpdateSingularByCross(*this,true);
            return true;
        }
        return false;
    }

    bool LoadSharpFeatures(const std::string &filename)
    {
        std::cout<<"Loading Sharp Features"<<std::endl;
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                face[i].FKind[j]=ETNone;
                face[i].ClearFaceEdgeS(j);
            }
        }

        FILE *f=fopen(filename.c_str(),"rt");
        if (f==NULL)return false;
        int Num;
        fscanf(f,"%d/n",&Num);

        for (size_t i=0;i<Num;i++)
        {
            int TypeSh,IndexF,IndexE;
            fscanf(f,"%d,%d,%d,/n",&TypeSh,&IndexF,&IndexE);
            assert((TypeSh==0)||(TypeSh==1));
            assert((IndexE>=0)&&(IndexE<3));
            assert((IndexF>=0)&&(IndexF<face.size()));
            if (TypeSh==0)
                face[IndexF].FKind[IndexE]=ETConcave;
            else
                face[IndexF].FKind[IndexE]=ETConvex;

            face[IndexF].SetFaceEdgeS(IndexE);

            if (!vcg::face::IsBorder(face[IndexF],IndexE))
            {
                FaceType *Fopp=face[IndexF].FFp(IndexE);
                int IOpp=face[IndexF].FFi(IndexE);
                Fopp->SetFaceEdgeS(IOpp);
                if (TypeSh==0)
                    Fopp->FKind[IOpp]=ETConcave;
                else
                    Fopp->FKind[IOpp]=ETConvex;
            }
        }
        fclose(f);
        return true;
    }

    bool LoadSharpFeaturesFL(const std::string &filename)
    {
        std::cout<<"Loading Sharp Features FL Format"<<std::endl;
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                face[i].FKind[j]=ETNone;
                face[i].ClearFaceEdgeS(j);
            }
        }

        FILE *f=fopen(filename.c_str(),"rt");
        if (f==NULL)return false;
        int Num;
        fscanf(f,"%d/n",&Num);

        std::map<std::pair<size_t,size_t> , std::pair<size_t,size_t> > FaceEdgeMap;
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                std::pair<size_t,size_t> FaceEdge(i,j);
                size_t IndexV0=vcg::tri::Index(*this,face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(*this,face[i].V1(j));
                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),
                                             std::max(IndexV0,IndexV1));
                FaceEdgeMap[key]=FaceEdge;
            }
        }

        for (size_t i=0;i<Num;i++)
        {
            int SizeSh;
            fscanf(f,"%d/n",&SizeSh);
            std::vector<int> CurrSh(SizeSh,-1);
            for (size_t j=0;j<SizeSh;j++)
                fscanf(f,"%d",&CurrSh[j]);

            for (size_t j=0;j<CurrSh.size()-1;j++)
            {
                size_t IndexV0=CurrSh[j];
                size_t IndexV1=CurrSh[j+1];
                //std::cout<<"Adding Sharp:"<<IndexV0<<";"<<IndexV1<<std::endl;
                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),
                                             std::max(IndexV0,IndexV1));
                assert(FaceEdgeMap.count(key)>0);
                size_t IndexF=FaceEdgeMap[key].first;
                size_t IndexE=FaceEdgeMap[key].second;
                face[IndexF].SetFaceEdgeS(IndexE);

                if (IsConcaveEdge(face[IndexF],IndexE))
                    face[IndexF].FKind[IndexE]=ETConcave;
                else
                    face[IndexF].FKind[IndexE]=ETConvex;

                FeatureKind  FKind=face[IndexF].FKind[IndexE];

                if (!vcg::face::IsBorder(face[IndexF],IndexE))
                {
                    FaceType *Fopp=face[IndexF].FFp(IndexE);
                    int IOpp=face[IndexF].FFi(IndexE);
                    Fopp->SetFaceEdgeS(IOpp);
                    Fopp->FKind[IOpp]=FKind;
                }
            }
        }
        fclose(f);
        return true;
    }

public:

    bool SaveTriMesh(const std::string &filename)
    {
        if(filename.empty()) return false;
        int position0=filename.find(".ply");
        int position1=filename.find(".obj");
        int position2=filename.find(".off");

        if (position0!=-1)
        {
            int err=vcg::tri::io::ExporterPLY<FieldTriMesh>::Save(*this,filename.c_str());
            if (err!=vcg::ply::E_NOERROR)return false;
            return true;
        }
        if (position1!=-1)
        {
            int mask=0;
            int err=vcg::tri::io::ExporterOBJ<FieldTriMesh>::Save(*this,filename.c_str(),mask);
            if ((err!=0)&&(err!=5))return false;
            return true;
        }
        if (position2!=-1)
        {
            int err=vcg::tri::io::ExporterOFF<FieldTriMesh>::Save(*this,filename.c_str());
            if (err!=0)return false;
            return true;
        }
        return false;
    }

    bool SaveOrigFace(const std::string &filename)
    {
        if(filename.empty()) return false;
        std::ofstream myfile;
        myfile.open (filename.c_str());

        for (size_t i=0;i<face.size();i++)
            myfile <<face[i].IndexOriginal <<std::endl;

        myfile.close();
        return true;
    }

    bool SaveField(const std::string &filename)
    {
        if(filename.empty()) return false;
        vcg::tri::io::ExporterFIELD<FieldTriMesh>::Save4ROSY(*this,filename.c_str());
        return true;
    }

    //VCG UPDATING STRUCTURES
    void UpdateDataStructures()
    {
        vcg::tri::Clean<FieldTriMesh>::RemoveUnreferencedVertex(*this);
        vcg::tri::Allocator<FieldTriMesh>::CompactEveryVector(*this);

        vcg::tri::UpdateBounding<FieldTriMesh>::Box(*this);
        vcg::tri::UpdateNormal<FieldTriMesh>::PerVertexNormalizedPerFace(*this);
        vcg::tri::UpdateNormal<FieldTriMesh>::PerFaceNormalized(*this);
        vcg::tri::UpdateTopology<FieldTriMesh>::FaceFace(*this);
        vcg::tri::UpdateTopology<FieldTriMesh>::VertexFace(*this);
        vcg::tri::UpdateFlags<FieldTriMesh>::VertexBorderFromFaceAdj(*this);
        vcg::tri::UpdateFlags<FieldTriMesh>::FaceBorderFromFF(*this);
    }


    void InitSharpFeatures(ScalarType SharpAngleDegree)
    {
        UpdateDataStructures();

        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
                face[i].ClearFaceEdgeS(j);

        if (SharpAngleDegree>0)
            vcg::tri::UpdateFlags<MeshType>::FaceEdgeSelCrease(*this,vcg::math::ToRad(SharpAngleDegree));
        InitEdgeType();

        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (vcg::face::IsBorder(face[i],j))
                {
                    face[i].SetFaceEdgeS(j);
                    face[i].FKind[j]=ETConvex;
                }
            }
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {

                if (!vcg::face::IsManifold(face[i],j))
                {
                    face[i].SetFaceEdgeS(j);
                    face[i].FKind[j]=ETConvex;
                }
            }
        std::cout<<"There is "<<SharpLenght()<<" sharp lenght"<<std::endl;
    }

    ScalarType SharpLenght()
    {
        ScalarType LSharp=0;
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                ScalarType L=(face[i].P0(j)-face[i].P1(j)).Norm();
                if (vcg::face::IsBorder(face[i],j))
                    LSharp+=2*L;
                else
                {
                    if (face[i].IsFaceEdgeS(j))
                        LSharp+=L;
                }
            }
        return (LSharp/2);
    }

    ScalarType Area()const
    {
        ScalarType CurrA=0;
        for (size_t i=0;i<face.size();i++)
            CurrA+=vcg::DoubleArea(face[i]);
        return (CurrA/2);
    }

    ScalarType SignedVolume(const size_t &IndexF)const
    {
        return (face[IndexF].cP(0)*(face[IndexF].cP(1)^face[IndexF].cP(2))/6.0);
    }

    ScalarType Volume()const
    {
        ScalarType Vol;
        for (size_t i=0;i<face.size();i++)
            Vol+=SignedVolume(i);

        return (fabs(Vol));
    }

    void SetFeatureValence()
    {
        vcg::tri::UpdateQuality<FieldTriMesh>::VertexConstant(*this,0);
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                face[i].V0(j)->Q()+=1;
                face[i].V1(j)->Q()+=1;
            }

    }

    void ErodeFeaturesStep()
    {
        SetFeatureValence();

        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                if (vcg::face::IsBorder(face[i],j))continue;
                ScalarType Len=(face[i].P0(j)-face[i].P1(j)).Norm();
                if (Len>bbox.Diag()*0.05)continue;

                if ((face[i].V0(j)->Q()==2)||(face[i].V1(j)->Q()==2))
                    face[i].ClearFaceEdgeS(j);
            }
    }

    void DilateFeaturesStep(std::vector<std::pair<size_t,size_t> > &OrigFeatures)
    {
        SetFeatureValence();

        for (size_t i=0;i<OrigFeatures.size();i++)
        {
            size_t IndexF=OrigFeatures[i].first;
            size_t IndexE=OrigFeatures[i].second;
            if ((face[IndexF].V0(IndexE)->Q()==2)&&
                    (!face[IndexF].V0(IndexE)->IsS()))
                face[IndexF].SetFaceEdgeS(IndexE);

            if ((face[IndexF].V1(IndexE)->Q()==2)&&
                    (!face[IndexF].V1(IndexE)->IsS()))
                face[IndexF].SetFaceEdgeS(IndexE);
        }
    }

    void PrintSharpInfo()
    {
        vcg::tri::UpdateQuality<MeshType>::VertexConstant(*this,0);
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                if (!face[i].V0(j)->IsB())
                    face[i].V0(j)->Q()+=1;
                if (!face[i].V1(j)->IsB())
                    face[i].V1(j)->Q()+=1;
            }
        size_t Valence1=0;
        for (size_t i=0;i<vert.size();i++)
            if (vert[i].Q()==2)Valence1++;
        std::cout<<"EndPoints "<<Valence1<<std::endl;
    }

    void ErodeDilate(size_t StepNum)
    {
        vcg::tri::UpdateFlags<FieldTriMesh>::VertexClearS(*this);
        SetFeatureValence();
        for (size_t i=0;i<vert.size();i++)
            if ((vert[i].Q()>4)||((vert[i].IsB())&&(vert[i].Q()>2)))
                vert[i].SetS();

        std::vector<std::pair<size_t,size_t> > OrigFeatures;

        //save the features
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                OrigFeatures.push_back(std::pair<size_t,size_t>(i,j));
            }

        for (size_t s=0;s<StepNum;s++)
            ErodeFeaturesStep();
        for (size_t s=0;s<StepNum;s++)
            DilateFeaturesStep(OrigFeatures);

        PrintSharpInfo();
    }

    bool SufficientFeatures(ScalarType SharpFactor)
    {
        ScalarType sqrtCurrA=sqrt(Area());
        ScalarType SharpL=SharpLenght();
        std::cout<<"Sqrt Area "<<sqrtCurrA<<std::endl;
        std::cout<<"Sharp Lenght "<<SharpL<<std::endl;
        ScalarType Ratio=SharpL/sqrtCurrA;
        std::cout<<"Ratio "<<Ratio<<std::endl;
        return(Ratio>SharpFactor);
    }

    void SetSharp(FaceType &f,int IndexE)
    {
        f.SetFaceEdgeS(IndexE);
        if (IsConcaveEdge(f,IndexE))
            f.FKind[IndexE]=ETConcave;
        else
            f.FKind[IndexE]=ETConvex;

        if (vcg::face::IsBorder(f,IndexE))return;

        FieldTriMesh::FaceType *fopp=f.FFp(IndexE);
        int eOpp=f.FFi(IndexE);

        fopp->SetFaceEdgeS(eOpp);
        fopp->FKind[eOpp]= f.FKind[IndexE];
    }

    void ClearSharp(FaceType &f,int IndexE)
    {
        f.ClearFaceEdgeS(IndexE);
        f.FKind[IndexE]=ETNone;
        if (vcg::face::IsBorder(f,IndexE))return;

        FieldTriMesh::FaceType *fopp=f.FFp(IndexE);
        int eOpp=f.FFi(IndexE);

        fopp->ClearFaceEdgeS(eOpp);
        fopp->FKind[eOpp]= ETNone;
    }

    void InitFaceOriginalIndex()
    {
        for (size_t i=0;i<face.size();i++)
            face[i].IndexOriginal=i;
    }
};


#endif
