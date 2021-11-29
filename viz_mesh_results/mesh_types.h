#ifndef DEFAULTMESHTYPES_H
#define DEFAULTMESHTYPES_H

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/import_field.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/color.h>

/* ----- Triangle mesh ----- */

class TriangleVertex;
class TriangleFace;
struct MyTriangleTypes : public vcg::UsedTypes<
        vcg::Use<TriangleVertex>::AsVertexType,
        vcg::Use<TriangleFace>::AsFaceType>{};

class TriangleVertex : public vcg::Vertex<
        MyTriangleTypes,
        vcg::vertex::Coord3d,
        vcg::vertex::Normal3d,
        vcg::vertex::VFAdj,
        vcg::vertex::Color4b,
        vcg::vertex::Qualityd,
        vcg::vertex::BitFlags,
        vcg::vertex::CurvatureDird>{};

class TriangleFace : public vcg::Face<
        MyTriangleTypes,
        vcg::face::VertexRef,
        vcg::face::Normal3d,
        vcg::face::Color4b,
        vcg::face::Qualityd,
        vcg::face::BitFlags,
        vcg::face::FFAdj,
        vcg::face::VFAdj,
        vcg::face::CurvatureDird,
        vcg::face::Mark,
        vcg::face::WedgeTexCoord2d> {};

class TriangleMesh : public vcg::tri::TriMesh<
        std::vector<TriangleVertex>,
        std::vector<TriangleFace> >
{
public:

    bool LoadMesh(std::string filename)
    {
        Clear();
        if(filename.empty()) return false;
        int position0=filename.find(".ply");
        int position1=filename.find(".obj");
        int position2=filename.find(".off");

        if (position0!=-1)
        {
            int err=vcg::tri::io::ImporterPLY<TriangleMesh>::Open(*this,filename.c_str());
            if (err!=vcg::ply::E_NOERROR)return false;
            return true;
        }
        if (position1!=-1)
        {
            int mask;
            vcg::tri::io::ImporterOBJ<TriangleMesh>::LoadMask(filename.c_str(),mask);
            int err=vcg::tri::io::ImporterOBJ<TriangleMesh>::Open(*this,filename.c_str(),mask);
            if ((err!=0)&&(err!=5))return false;
            return true;
        }
        if (position2!=-1)
        {
            int err=vcg::tri::io::ImporterOFF<TriangleMesh>::Open(*this,filename.c_str());
            if (err!=0)return false;
            return true;
        }
        return false;
    }

    bool LoadField(std::string field_filename)
    {
        int position0=field_filename.find(".ffield");
        int position1=field_filename.find(".rosy");


        if (position0!=-1)
        {
            bool loaded=vcg::tri::io::ImporterFIELD<TriangleMesh>::LoadFFIELD(*this,field_filename.c_str());
            if (!loaded)return false;
            vcg::tri::CrossField<TriangleMesh>::OrientDirectionFaceCoherently(*this);
            vcg::tri::CrossField<TriangleMesh>::UpdateSingularByCross(*this,true);
            return true;
        }
        if (position1!=-1)
        {
            std::cout<<"Importing ROSY field"<<std::endl;
            bool loaded=vcg::tri::io::ImporterFIELD<TriangleMesh>::Load4ROSY(*this,field_filename.c_str());
            std::cout<<"Imported ROSY field"<<std::endl;
            if (!loaded)return false;
            vcg::tri::CrossField<TriangleMesh>::OrientDirectionFaceCoherently(*this);
            vcg::tri::CrossField<TriangleMesh>::UpdateSingularByCross(*this,true);
            return true;
        }
        return false;
    }

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

    bool LoadSharpFeatures(std::string &FeaturePath)
    {
        std::cout<<"Loading Sharp Features"<<std::endl;
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
        }
        fclose(f);
        std::cout<<"Loaded Sharp Features"<<std::endl;
        return true;
    }



    void loadPatchesIntoQ(const std::string& filename)
    {
        std::ifstream input;
        input.open(filename.c_str());
        if (!input.is_open())
        {
            std::cout<<"ERROR LOADING PATCH FILE"<<std::endl;
            exit(0);
        }

        size_t numFaces;
        input >> numFaces;

        for (size_t i=0; i < numFaces;i++)
        {
            int CurrP;
            input >> CurrP;
            face[i].Q()=CurrP;
        }

        input.close();
    }

    void GenerateSingMesh(TriangleMesh &SingMesh,
                          ScalarType scale=0.005)
    {
       SingMesh.Clear();
        // query if an attribute is present or not
       bool hasSingular = vcg::tri::HasPerVertexAttribute(*this,std::string("Singular"));
       bool hasSingularIndex = vcg::tri::HasPerVertexAttribute(*this,std::string("SingularIndex"));

       if (!hasSingular)return;
       if(!hasSingularIndex)return;

       typename MeshType::template PerVertexAttributeHandle<bool> Handle_Singular;
       Handle_Singular=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<bool>(*this,std::string("Singular"));
       typename MeshType::template PerVertexAttributeHandle<int> Handle_SingularIndex;
       Handle_SingularIndex =vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(*this,std::string("SingularIndex"));

       ScalarType sizeP=bbox.Diag()*scale;
       for (size_t i=0;i<vert.size();i++)
       {
           if (vert[i].IsD())continue;
           if (!Handle_Singular[i])continue;


           int SingIndex=Handle_SingularIndex[i];

           vcg::Color4b colSing;

           switch (SingIndex)
           {
             case 1:colSing=vcg::Color4b(0,0,255,255);      break;
             case 2:colSing=vcg::Color4b(0,255,0,255);    break;
             case 3:colSing=vcg::Color4b(255,0,0,255);      break;
             case 4:colSing=vcg::Color4b(255,255,0,255);      break;
             default:colSing=vcg::Color4b(255,0,255,255);
           }

           TriangleMesh vertMesh;
           vcg::tri::Sphere<TriangleMesh>(vertMesh);
           for (size_t j=0;j<vertMesh.vert.size();j++)
           {
               vertMesh.vert[j].P()*=sizeP;
               vertMesh.vert[j].P()+=vert[i].P();
           }

           vcg::tri::UpdateColor<TriangleMesh>::PerFaceConstant(vertMesh,colSing);
           vcg::tri::Append<TriangleMesh,TriangleMesh>::Mesh(SingMesh,vertMesh);
       }
    }

    void GenerateFieldMesh(TriangleMesh &FieldMesh,
                           ScalarType scale=0.01)
    {
        FieldMesh.Clear();
        ScalarType sizeF=bbox.Diag()*scale;
        for (size_t i=0;i<face.size();i++)
        {
            CoordType Pos0=(face[i].P(0)+face[i].P(1)+face[i].P(2))/3;
            CoordType Dir[4];
            Dir[0]=face[i].PD1();
            Dir[1]=face[i].PD2();
            Dir[2]=-face[i].PD1();
            Dir[3]=-face[i].PD2();
            CoordType Pos[4];
            for (size_t j=0;j<4;j++)
                Pos[j]=Pos0+Dir[j]*sizeF;

            for (size_t j=0;j<4;j++)
            {
                TriangleMesh cylMesh;
                vcg::tri::OrientedCylinder<TriangleMesh>(cylMesh, Pos0, Pos[j], sizeF/20, true,4,4);
                vcg::tri::UpdateColor<TriangleMesh>::PerFaceConstant(cylMesh,vcg::Color4b::Black);
                vcg::tri::Append<TriangleMesh,TriangleMesh>::Mesh(FieldMesh,cylMesh);
            }
        }
        FieldMesh.UpdateAttributes();
    }

    void GenerateEdgeSelMesh(TriangleMesh &SharpMesh,ScalarType scale=0.01)
    {
        SharpMesh.Clear();
        ScalarType sizeS=bbox.Diag()*scale;
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                bool IsB=vcg::face::IsBorder(face[i],j);
                bool IsS=face[i].IsFaceEdgeS(j);
                if (!IsB)
                    IsS|=face[i].FFp(j)->IsFaceEdgeS(face[i].FFi(j));

                bool Order=(face[i].V0(j)<face[i].V1(j));
                bool AddCyl=IsB;
                AddCyl|=(IsS && Order);
                if (!AddCyl)continue;
                CoordType Pos0=face[i].P0(j);
                CoordType Pos1=face[i].P1(j);
                TriangleMesh cylMesh;
                vcg::tri::OrientedCylinder<TriangleMesh>(cylMesh, Pos0, Pos1, sizeS/10, true,4,4);
                vcg::tri::UpdateColor<TriangleMesh>::PerFaceConstant(cylMesh,vcg::Color4b::Green);
                vcg::tri::Append<TriangleMesh,TriangleMesh>::Mesh(SharpMesh,cylMesh);
            }
        }
        SharpMesh.UpdateAttributes();
    }

    void SelectPatchBorders()
    {
        vcg::tri::UpdateSelection<TriangleMesh>::FaceClear(*this);
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                bool IsB=vcg::face::IsBorder(face[i],j);
                if (IsB)
                {
                    face[i].SetFaceEdgeS(j);
                    continue;
                }
                int Q0=face[i].Q();
                int Q1=face[i].FFp(j)->Q();
                if (Q0!=Q1)
                    face[i].SetFaceEdgeS(j);
            }
        }
    }

    void ColorByPartition()
    {
        int MaxQ=0;
        for (size_t i=0;i<face.size();i++)
            MaxQ=std::max(MaxQ,(int)face[i].Q());
        for (size_t i=0;i<face.size();i++)
            face[i].C()=vcg::Color4b::Scatter(MaxQ,face[i].Q());
    }
};


#endif // DEFAULTMESHTYPES_H
