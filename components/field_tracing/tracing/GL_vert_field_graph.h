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

#ifndef GL_VERT_FIELD_GRAPH
#define GL_VERT_FIELD_GRAPH

#include "vert_field_graph.h"

template <class MeshType>
class GLVertGraph
{

    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;

    VertexFieldGraph<MeshType> &VFGraph;

    std::vector<CoordType> DisplacedP;


public:

    void GLDrawPoints(const std::vector<CoordType> &DrawPos,
                      const ScalarType &GLSize,const vcg::Color4b &Col)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9995);
        glPointSize(GLSize);

        for (size_t i=0;i<DrawPos.size();i++)
        {

            vcg::glColor(Col);
            glBegin(GL_POINTS);
            vcg::glVertex(DrawPos[i]);
            glEnd();
        }
        glPopAttrib();
    }

    void GLDrawAllNodesNeigh()
    {
        for (size_t i=0;i<VFGraph.NumNodes();i++)
            GLDrawNodeNeigh(i,5);
    }

    void GLDrawNodes(const std::vector<size_t> &IndexN,
                     const ScalarType &size,
                     bool DrawDir=true,
                     ScalarType GLSize=10,
                     bool DrawNeigh=false)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9995);
        glPointSize(GLSize);
        glLineWidth(GLSize/2);

        for (size_t i=0;i<IndexN.size();i++)
        {
            size_t currN=IndexN[i];

            CoordType Dir=VFGraph.NodeDir(currN);

            vcg::glColor(vcg::Color4b(0,0,255,255));
            glBegin(GL_POINTS);
            vcg::glVertex(DisplacedP[currN]);
            glEnd();

            vcg::glColor(vcg::Color4b(0,255,0,255));
            if (DrawDir)
            {
                glBegin(GL_LINES);
                vcg::glVertex(DisplacedP[currN]);
                vcg::glVertex(DisplacedP[currN]+Dir*size/2);
                glEnd();
            }
            if (DrawNeigh)
                GLDrawNodeNeigh(currN,2);
        }

        glPopAttrib();
    }

    void GLDrawNonActiveNodes(const ScalarType &size,bool DrawDir=true)
    {
        std::vector<size_t> ActiveList;
        for (size_t i=0;i<VFGraph.NumNodes();i++)
            if (!VFGraph.IsActive(i))ActiveList.push_back(i);
        GLDrawNodes(ActiveList,size,DrawDir);
    }

    void GLDrawSingNodes(const ScalarType &size,
                         bool DrawDir=true,
                         ScalarType GLSize=10)
    {
        GLDrawNodes(VFGraph.SingNodes,size,DrawDir,GLSize);
    }

    //    void GLDrawNode(size_t &IndexNode,
    //                    const vcg::Color4b &NodeCol,
    //                    const ScalarType &PointSize=10,
    //                    const ScalarType &Size=0,
    //                    const CoordType &Displ=CoordType(0,0,0),
    //                    const bool MoveOnDir=false)
    //    {
    //        CoordType Pos,Dir;
    //        VFGraph.NodePosDir(IndexNode,Pos,Dir);
    //        vcg::glColor(NodeCol);

    //        if (MoveOnDir)
    //            Pos+=Dir*Size/2;

    //        glPointSize(PointSize);

    //        glPushAttrib(GL_ALL_ATTRIB_BITS);
    //        glDisable(GL_LIGHTING);
    //        glDepthRange(0,0.9998);

    //        glBegin(GL_POINTS);
    //        vcg::glVertex(Pos+Displ);
    //        glEnd();

    //        if (Size>0)
    //        {
    //            glLineWidth(std::max(PointSize/5.0,1.0));
    //            glBegin(GL_LINES);
    //            vcg::glVertex(Pos);
    //            vcg::glVertex(Pos+Dir*Size);
    //            glEnd();
    //        }

    //        glPopAttrib();
    //    }

    void InitDisplacedPos()
    {
        //get all average positions
        DisplacedP.clear();
        std::vector<CoordType> AvgPos(VFGraph.Mesh().vert.size(),CoordType(0,0,0));
        std::vector<size_t> AvgNum(VFGraph.Mesh().vert.size(),0);
        for (size_t i=0;i<VFGraph.Mesh().face.size();i++)
        {
            CoordType Bary=(VFGraph.Mesh().face[i].P(0)+
                            VFGraph.Mesh().face[i].P(1)+
                            VFGraph.Mesh().face[i].P(2))/3;
            for (size_t j=0;j<3;j++)
            {
                VertexType *v=VFGraph.Mesh().face[i].V0(j);
                size_t IndexV=vcg::tri::Index(VFGraph.Mesh(),v);
                AvgPos[IndexV]+=Bary;
                AvgNum[IndexV]++;
            }
        }

        std::vector<CoordType> TargetP;
        for (size_t i=0;i<AvgPos.size();i++)
        {
            AvgPos[i]/=AvgNum[i];
            CoordType PosReal=VFGraph.Mesh().vert[i].P();
            for (size_t j=0;j<4;j++)
                DisplacedP.push_back(PosReal*0.5+AvgPos[i]*0.5);
        }

        ScalarType size=VFGraph.Mesh().bbox.Diag()*0.0005;

        for (size_t i=0;i<DisplacedP.size();i++)
        {
            CoordType Dir=VFGraph.NodeDir(i);
            DisplacedP[i]+=Dir*size/2;
        }
    }

    void GLDrawNodeNeigh(const size_t IndexNode,
                         const ScalarType GLSize)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9995);
        glLineWidth(GLSize);

        CoordType Pos0=DisplacedP[IndexNode];
        for (size_t i=0;i<VFGraph.NumNeigh(IndexNode);i++)
        {
            size_t IndexNeigh=VFGraph.NodeNeigh(IndexNode,i);
            CoordType Pos1=DisplacedP[IndexNeigh];
            if (VFGraph.DirectNeigh(IndexNode,i))
                vcg::glColor(vcg::Color4b(0,255,255,255));
            else
                continue;

            glBegin(GL_LINES);
            vcg::glVertex(Pos0);
            vcg::glVertex(Pos1);
            glEnd();
        }

        glPopAttrib();
    }

    void GLDrawTwinsConnections()
    {

        ScalarType size=VFGraph.Mesh().bbox.Diag()*0.0005;

        std::vector<size_t> IndexN;
        for (size_t i=0;i<VFGraph.NumNodes();i++)
        {
            if (!VFGraph.HasTwin(i))continue;
            IndexN.push_back(i);
        }
        GLDrawNodes(IndexN,size);
        //CoordType Dir=VFGraph.NodeDir(i);

        //            vcg::glColor(vcg::Color4b(0,0,255,255));
        //            glBegin(GL_POINTS);
        //            vcg::glVertex(DisplacedP[i]);
        //            glEnd();

        //            vcg::glColor(vcg::Color4b(0,255,0,255));
        ////            if (DrawDir)
        ////            {
        //                glBegin(GL_LINES);
        //                vcg::glVertex(DisplacedP[i]);
        //                vcg::glVertex(DisplacedP[i]+Dir*size/2);
        //                glEnd();
        //            //}

        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.995);
        glPointSize(10);
        glLineWidth(5);

        for (size_t i=0;i<VFGraph.NumNodes();i++)
        {
            if (!VFGraph.HasTwin(i))continue;
            for (size_t j=0;j<VFGraph.NumNeigh(i);j++)
            {
                if (!VFGraph.TwinNeigh(i,j))continue;
                size_t IndexN0=i;
                size_t IndexN1=VFGraph.NodeNeigh(i,j);
                CoordType Pos0=DisplacedP[IndexN0];
                CoordType Pos1=DisplacedP[IndexN1];
                glBegin(GL_LINES);
                vcg::glColor(vcg::Color4b(255,255,0,255));
                vcg::glVertex(Pos0);
                vcg::glColor(vcg::Color4b(255,0,0,255));
                vcg::glVertex(Pos1);
                glEnd();
            }
        }

        glPopAttrib();
    }

    void GLDrawPath(const std::vector<size_t> &Path,
                    vcg::Color4b PathCol,
                    bool IsLoop)
    {
        vcg::glColor(PathCol);

        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9998);
        glLineWidth(12);
        glBegin(GL_LINE_STRIP);
        for (size_t i=0;i<Path.size();++i)
        {
            CoordType Pos=VFGraph.NodePos(Path[i]);
            vcg::glVertex(Pos);
        }
        if (IsLoop)
            vcg::glVertex(VFGraph.NodePos(Path[0]));


        glEnd();

        glPopAttrib();
    }

    void GLDrawPaths(const std::vector<std::vector<size_t> > &Paths,
                     const std::vector<bool> &IsLoop,
                     ScalarType size,
                     bool DrawNode,
                     int VisTraces=-1,
                     bool ColorByLoop=false)
    {
        size_t Limit=Paths.size();
        if ((VisTraces>0)&&(VisTraces<Paths.size()))
            Limit=VisTraces;
        for (size_t i=0;i<Limit;i++)
        {
            vcg::Color4b Col=vcg::Color4b::Scatter(Paths.size(),i);

            if (ColorByLoop)
            {
                if (IsLoop[i])
                    Col=vcg::Color4b::Blue;
                else
                    Col=vcg::Color4b::Red;
            }

            GLDrawPath(Paths[i],Col,IsLoop[i]);
            if (DrawNode)
                GLDrawNodes(Paths[i],size);
            //            for (size_t j=0;j<Paths[i].size();j++)
            //                GLDrawNodeNeigh(Paths[i][j],10);
        }
    }

    GLVertGraph(VertexFieldGraph<MeshType> &_VFGraph):VFGraph(_VFGraph)
    {

    }
};

#endif
