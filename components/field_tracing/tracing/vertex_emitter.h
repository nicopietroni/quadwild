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

#ifndef VERTEX_EMITTER
#define VERTEX_EMITTER

#include "vert_field_graph.h"

template <class MeshType>
class VertexEmitter
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;


public:

    static void ComputeFlatEmitterReceivers(VertexFieldGraph<MeshType> &VFGraph,
                                            const std::vector<std::vector<CoordType> > &VertOrthoDir,
                                            //const std::vector<std::vector<CoordType> > &VertFlatDir,
                                            const size_t IndexV,
                                            size_t &Emitter,
                                            size_t &Receiver)
    {
        //assert(VertOrthoDir[IndexV].size()==2);
        if (VertOrthoDir[IndexV].size()!=2)
        {
           std::cout<<"WARNING ON BORDER There are "<<VertOrthoDir[IndexV].size()<<" Ortho dirs instead of 2"<<std::endl;
           //assert(0);
        }
        CoordType Ortho0=VertOrthoDir[IndexV][0];
        CoordType Ortho1=VertOrthoDir[IndexV][1];
        CoordType TargetD=Ortho0+Ortho1;
        TargetD.Normalize();
        size_t BestDir=VFGraph.GetClosestDirTo(IndexV,TargetD);
        Emitter=VertexFieldGraph<MeshType>::IndexNode(IndexV,BestDir);
        Receiver=VertexFieldGraph<MeshType>::TangentNode(Emitter);
    }

    static void ComputeNarrowEmitterReceivers(VertexFieldGraph<MeshType> &VFGraph,
                                              //const std::vector<std::vector<CoordType> > &VertOrthoDir,
                                              const std::vector<std::vector<CoordType> > &VertFlatDir,
                                              const size_t IndexV,
                                              size_t &Emitter,
                                              size_t &Receiver)
    {
        assert(VertFlatDir[IndexV].size()==2);
        CoordType Flat0=VertFlatDir[IndexV][0];
        CoordType Flat1=VertFlatDir[IndexV][1];
        CoordType TargetD=Flat0+Flat1;
        TargetD.Normalize();
        size_t BestDir=VFGraph.GetClosestDirTo(IndexV,TargetD);
        Emitter=VertexFieldGraph<MeshType>::IndexNode(IndexV,BestDir);
        Receiver=VertexFieldGraph<MeshType>::TangentNode(Emitter);
    }

    static void ComputeConcaveEmitterReceivers(VertexFieldGraph<MeshType> &VFGraph,
                                               const std::vector<std::vector<CoordType> > &VertOrthoDir,
                                               const std::vector<std::vector<CoordType> > &VertFlatDir,
                                               const size_t IndexV,
                                               std::vector<size_t> &Emitter,
                                               std::vector<size_t> &Receiver)
    {
        assert(VertOrthoDir[IndexV].size()==2);
//        if(VertOrthoDir[IndexV].size()!=2)
//        {
//            std::cout<<VertOrthoDir[IndexV].size()<<std::endl;
//            for (size_t j=0;j<VertOrthoDir[IndexV].size();j++)
//            {
//                CoordType Dir=VertOrthoDir[IndexV][j];
//                std::cout<<"Dir:"<<Dir.X()<<","<<Dir.Y()<<","<<Dir.Z()<<std::endl;
//            }
//            assert(0);
//        }
        CoordType Ortho0=VertOrthoDir[IndexV][0];
        CoordType Ortho1=VertOrthoDir[IndexV][1];
        size_t BestDir0=VFGraph.GetClosestDirTo(IndexV,Ortho0);
        size_t BestDir1=VFGraph.GetClosestDirTo(IndexV,Ortho1);

        size_t IndexNode0=VertexFieldGraph<MeshType>::IndexNode(IndexV,BestDir0);
        size_t IndexNode1=VertexFieldGraph<MeshType>::IndexNode(IndexV,BestDir1);
        //there are three cases
        if (IndexNode0==IndexNode1)
            Emitter.push_back(IndexNode0);
        else
        {
            Emitter.push_back(IndexNode0);
            Emitter.push_back(IndexNode1);
            //check if they are opposite
            size_t OppN=VertexFieldGraph<MeshType>::TangentNode(IndexNode0);
            if (OppN==IndexNode1)//in this case add a third one
            {
                size_t EmitterN,ReceiverN;
                //ComputeNarrowEmitterReceivers(VFGraph,VertOrthoDir,VertFlatDir,IndexV,EmitterN,ReceiverN);
                ComputeNarrowEmitterReceivers(VFGraph,VertFlatDir,IndexV,EmitterN,ReceiverN);
                if ((EmitterN!=IndexNode0)&&(EmitterN!=IndexNode1))
                    Emitter.push_back(EmitterN);
                else
                {
                    std::cout<<"WARNING THIRD emitter problem"<<std::endl;
//                    vcg::tri::UpdateFlags<MeshType>::VertexClear(VFGraph.Mesh());
//                    VFGraph.Mesh().vert[IndexV].SetS();
//                    vcg::tri::io::ExporterPLY<MeshType>::Save(VFGraph.Mesh(),"test1.ply",vcg::tri::io::Mask::IOM_FLAGS);
//                    assert(0);
//                    assert(EmitterN!=IndexNode0);
//                    assert(EmitterN!=IndexNode1);
                }
            }
        }

        for (size_t i=0;i<Emitter.size();i++)
            Receiver.push_back(VertexFieldGraph<MeshType>::TangentNode(Emitter[i]));
    }

//    static void GetOrthoFlatDirections(VertexFieldGraph<MeshType> &VFGraph,
//                                       std::vector<std::vector<CoordType> > &VertFlatDir,
//                                       std::vector<std::vector<CoordType> > &VertOrthoDir)
//    {
//        VertOrthoDir.clear();
//        VertOrthoDir.resize(VFGraph.Mesh().vert.size());
//        VertFlatDir.clear();
//        VertFlatDir.resize(VFGraph.Mesh().vert.size());
//        for (size_t i=0;i<VFGraph.Mesh().face.size();i++)
//            for (size_t j=0;j<3;j++)
//            {
//                if(!vcg::face::IsBorder(VFGraph.Mesh().face[i],j))continue;

//                size_t IndexV0=vcg::tri::Index(VFGraph.Mesh(),VFGraph.Mesh().face[i].cV0(j));
//                size_t IndexV1=vcg::tri::Index(VFGraph.Mesh(),VFGraph.Mesh().face[i].cV1(j));
//                CoordType P0=VFGraph.Mesh().vert[IndexV0].P();
//                CoordType P1=VFGraph.Mesh().vert[IndexV1].P();
//                CoordType Dir=P1-P0;
//                Dir.Normalize();

//                CoordType OrthoDir=Dir;
//                OrthoDir=VFGraph.Mesh().face[i].cN()^OrthoDir;

//                vcg::Matrix33<ScalarType> Rot0=vcg::RotationMatrix(VFGraph.Mesh().face[i].cN(),VFGraph.Mesh().vert[IndexV0].N());
//                vcg::Matrix33<ScalarType> Rot1=vcg::RotationMatrix(VFGraph.Mesh().face[i].cN(),VFGraph.Mesh().vert[IndexV1].N());

//                CoordType OrthoDir0=Rot0*OrthoDir;
//                CoordType OrthoDir1=Rot1*OrthoDir;

//                OrthoDir0.Normalize();
//                OrthoDir1.Normalize();

//                VertOrthoDir[IndexV0].push_back(OrthoDir0);
//                VertOrthoDir[IndexV1].push_back(OrthoDir1);
//                VertFlatDir[IndexV0].push_back(-Dir);
//                VertFlatDir[IndexV1].push_back(Dir);
//            }
//    }

    static void GetOrthoFlatDirections(const MeshType &mesh,
                                       std::vector<std::vector<CoordType> > &VertFlatDir,
                                       std::vector<std::vector<CoordType> > &VertOrthoDir)
    {
        VertOrthoDir.clear();
        VertOrthoDir.resize(mesh.vert.size());
        VertFlatDir.clear();
        VertFlatDir.resize(mesh.vert.size());
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if(!vcg::face::IsBorder(mesh.face[i],j))continue;

                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].cV0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].cV1(j));
                CoordType P0=mesh.vert[IndexV0].P();
                CoordType P1=mesh.vert[IndexV1].P();
                CoordType Dir=P1-P0;
                Dir.Normalize();

                CoordType OrthoDir=Dir;
                OrthoDir=mesh.face[i].cN()^OrthoDir;

                vcg::Matrix33<ScalarType> Rot0=vcg::RotationMatrix(mesh.face[i].cN(),mesh.vert[IndexV0].N());
                vcg::Matrix33<ScalarType> Rot1=vcg::RotationMatrix(mesh.face[i].cN(),mesh.vert[IndexV1].N());

                CoordType OrthoDir0=Rot0*OrthoDir;
                CoordType OrthoDir1=Rot1*OrthoDir;

                OrthoDir0.Normalize();
                OrthoDir1.Normalize();

                VertOrthoDir[IndexV0].push_back(OrthoDir0);
                VertOrthoDir[IndexV1].push_back(OrthoDir1);
                VertFlatDir[IndexV0].push_back(-Dir);
                VertFlatDir[IndexV1].push_back(Dir);
            }
    }

};

#endif
