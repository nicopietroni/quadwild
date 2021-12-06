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

#ifndef MESH_FIELD_SMOOTHER
#define MESH_FIELD_SMOOTHER

#include <wrap/io_trimesh/export_field.h>
#include <wrap/io_trimesh/import.h>
#include <fields/field_smoother.h>

template <class MeshType>
class MeshFieldSmoother
{

    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType  CoordType;
    typedef typename MeshType::VertexType    VertexType;
    typedef typename MeshType::VertexPointer VertexPointer;
    typedef typename MeshType::FaceType    FaceType;
    typedef typename MeshType::FacePointer FacePointer;

    typedef vcg::tri::FieldSmoother<MeshType> FieldSmootherType;
    typedef typename FieldSmootherType::SmoothParam SmoothParam;

   public:

    static void AutoSetupParam(MeshType &mesh,typename vcg::tri::FieldSmoother<MeshType>::SmoothParam &FParam,ScalarType SharpFactor=6)
    {
        FParam.align_borders=true;
        FParam.sharp_thr=0;

        bool SufficientFeatures=mesh.SufficientFeatures(SharpFactor);

        if (SufficientFeatures)
        {
            std::cout<<"Using NPoly"<<std::endl;
            FParam.SmoothM=vcg::tri::SMNPoly;
            FParam.curv_thr=0;
            FParam.alpha_curv=0;
        }
        else
        {
#ifndef COMISO_FIELD
            std::cout<<"Using NPoly"<<std::endl;
            FParam.SmoothM=vcg::tri::SMNPoly;
            FParam.curv_thr=0.8;
            FParam.alpha_curv=0;
#else
            std::cout<<"Using Comiso"<<std::endl;
            FParam.SmoothM=vcg::tri::SMMiq;
            if (FParam.alpha_curv==0)
                FParam.alpha_curv=0.3;
            FParam.curv_thr=0;
#endif
        }
    }

    static void SmoothField(MeshType &mesh,SmoothParam FieldParam)
    {
        //CHANGE BUT NOT REALLY NEEDED
        //InitFeatureCoordsTable();

        FieldParam.sharp_thr=0.0;
        //FieldParam.curv_thr=0.0;
        FieldParam.AddConstr.clear();
        for (size_t i=0;i<mesh.face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                if (!mesh.face[i].IsFaceEdgeS(j))continue;
                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                //only on one side
                if (IndexV0>IndexV1)continue;
                CoordType P0=mesh.face[i].V0(j)->P();
                CoordType P1=mesh.face[i].V1(j)->P();
                CoordType Dir=P1-P0;
                Dir.Normalize();
                FieldParam.AddConstr.push_back(std::pair<int,CoordType>(i,Dir));

                typename MeshType::FaceType *FOpp=mesh.face[i].FFp(j);
                if (FOpp==&mesh.face[i])continue;

                int IndexF=vcg::tri::Index(mesh,*FOpp);
                FieldParam.AddConstr.push_back(std::pair<int,CoordType>(IndexF,Dir));
            }
        }
//        //if there is no alpha curvature and no constraint then set curv_thr
//        if ((FieldParam.alpha_curv==0)&&(FieldParam.AddConstr.size()==0))
//        {
//            FieldParam.curv_thr=0.1;
            FieldParam.align_borders=true;
//        }
#ifndef COMISO_FIELD
         FieldParam.SmoothM=vcg::tri::SMNPoly;
#endif
//        std::cout<<"..Alpha.."<<FieldParam.alpha_curv<<std::endl;
//        std::cout<<"..Smoothing.."<<std::endl;

//        vcg::tri::io::ExporterPLY<MeshType>::Save(mesh,"test.ply");
//        vcg::tri::io::ImporterPLY<MeshType>::Open(mesh,"test.ply");
//        mesh.UpdateDataStructures();

        FieldSmootherType::SmoothDirections(mesh,FieldParam);


        std::cout<<"..Global Orientation.."<<std::endl;
        vcg::tri::CrossField<MeshType>::OrientDirectionFaceCoherently(mesh);
        std::cout<<"..Update Singularities.."<<std::endl;
        vcg::tri::CrossField<MeshType>::UpdateSingularByCross(mesh);
        std::cout<<"..Update Features.."<<std::endl;


        //CHANGE BUT NOT REALLY NEEDED
        mesh.UpdateDataStructures();
        mesh.SetFeatureFromTable();
    }


    static bool SaveField(MeshType &mesh,const std::string &filename)
    {
        if(filename.empty()) return false;
        vcg::tri::io::ExporterFIELD<MeshType>::Save4ROSY(mesh,filename.c_str());
        return true;
    }

};


#endif
