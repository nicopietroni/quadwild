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

#ifndef MESH_DRAWING_H
#define MESH_DRAWING_H

#include <wrap/gl/trimesh.h>
#include <vcg/complex/algorithms/parametrization/uv_utils.h>
#include <vcg/complex/algorithms/update/topology.h>

template <class MeshType>
class MeshDrawing
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename vcg::Point2<ScalarType> UVCoordType;

    static CoordType UVTo3DPos(const vcg::Point2<ScalarType> &UVPos)
    {
        CoordType Pos(0,0,0);
        Pos.X()=UVPos.X();
        Pos.Y()=UVPos.Y();
        return Pos;
    }

public:

    static void GLDrawSharpEdges(MeshType &m,
                                 vcg::Color4b col=vcg::Color4b(255,0,0,255),
                                 ScalarType GLSize=5)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9999);
        glLineWidth(GLSize);
        vcg::glColor(col);
        glBegin(GL_LINES);
        for (size_t i=0;i<m.SharpFeatures.size();i++)
        {

            size_t IndexF=m.SharpFeatures[i].first;
            size_t IndexE=m.SharpFeatures[i].second;
            CoordType Pos0=m.face[IndexF].cP0(IndexE);
            CoordType Pos1=m.face[IndexF].cP1(IndexE);
            vcg::glVertex(Pos0);
            vcg::glVertex(Pos1);
        }
        glEnd();
        glPopAttrib();
    }


    static void GLDrawEdgeUV(MeshType &m)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glPushMatrix();
        vcg::Box2<ScalarType> uv_box=vcg::tri::UV_Utils<MeshType>::PerWedgeUVBox(m);
        vcg::glScale(3.0f/uv_box.Diag());
        vcg::glTranslate(CoordType(-uv_box.Center().X(),
                                   -uv_box.Center().Y(),0));
        glDisable(GL_LIGHTING);
        glDisable(GL_LIGHT0);
        glLineWidth(5);
        //glDepthRange(0,0.9999);
        glBegin(GL_LINES);
        for (size_t i=0;i<m.face.size();i++)
        {
            CoordType Pos[3];
            Pos[0]=UVTo3DPos(m.face[i].WT(0).P());
            Pos[1]=UVTo3DPos(m.face[i].WT(1).P());
            Pos[2]=UVTo3DPos(m.face[i].WT(2).P());
            for (size_t j=0;j<3;j++)
            {
                bool IsB=vcg::face::IsBorder(m.face[i],j);
                bool IsS=m.face[i].IsFaceEdgeS(j);
                if (!(IsS || IsB)) continue;
                vcg::glColor(vcg::Color4b(0,0,0,255));
                vcg::glVertex(Pos[j]);
                vcg::glVertex(Pos[(j+1)%3]);
            }
        }
        glEnd();
        glPopMatrix();
        glPopAttrib();
    }

    static void GLDrawUVPolylines(MeshType &m,
                                  std::vector<std::vector<UVCoordType> > &UVPolyL,
                                  std::vector<vcg::Color4b> &Color,
                                  std::vector<UVCoordType> &Dots,
                                  const vcg::Color4b &ColorDots)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glPushMatrix();
        vcg::Box2<ScalarType> uv_box=vcg::tri::UV_Utils<MeshType>::PerWedgeUVBox(m);
        vcg::glScale(3.0f/uv_box.Diag());
        //ScalarType UVScale=3.0f/uv_box.Diag();
        vcg::glTranslate(CoordType(-uv_box.Center().X(),
                                   -uv_box.Center().Y(),0));
        glDisable(GL_LIGHTING);
        glDisable(GL_LIGHT0);
        glLineWidth(10);


        for (size_t i=0;i<UVPolyL.size();i++)
        {
            if(UVPolyL[i].size()<2)continue;

            glBegin(GL_LINES);

            //std::cout<<"Size:"<<UVPolyL[i].size()<<std::endl;
            vcg::glColor(Color[i]);
            for (size_t j=0;j<UVPolyL[i].size()-1;j++)
            {
                CoordType P0(UVPolyL[i][j].X(),UVPolyL[i][j].Y(),0);
                CoordType P1(UVPolyL[i][j+1].X(),UVPolyL[i][j+1].Y(),0);
                vcg::glVertex(P0);
                vcg::glVertex(P1);
            }

            glEnd();
        }

        glPointSize(30);
        vcg::glColor(ColorDots);
        glBegin(GL_POINTS);
        for (size_t i=0;i<Dots.size();i++)
            vcg::glVertex(Dots[i]);
        glEnd();

        glPopMatrix();
        glPopAttrib();
    }

    static void GLDrawUV(MeshType &m,
                  int TxtIndex=-1,
                  bool colorPerVert=false)
    {
//        vcg::Box2<ScalarType> uv_box=vcg::tri::UV_Utils<CMesh>::PerWedgeUVBox(deformed_mesh);
//        ScalarType UVScale=3.0f/uv_box.Diag();
//        for (size_t i=0;i<deformed_mesh.face.size();i++)
//        {
//            deformed_mesh.face[i].WT(0).P()*=UVScale;
//            deformed_mesh.face[i].WT(1).P()*=UVScale;
//            deformed_mesh.face[i].WT(2).P()*=UVScale;
//        }

        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glPushMatrix();
        vcg::Box2<ScalarType> uv_box=vcg::tri::UV_Utils<MeshType>::PerWedgeUVBox(m);
        vcg::glScale(3.0f/uv_box.Diag());
        //ScalarType UVScale=3.0f/uv_box.Diag();
        vcg::glTranslate(CoordType(-uv_box.Center().X(),
                                   -uv_box.Center().Y(),0));
        glDisable(GL_LIGHTING);
        glDisable(GL_LIGHT0);

        if (TxtIndex>=0)
        {
            glActiveTexture(GL_TEXTURE0);
            glEnable(GL_TEXTURE_2D);
            glBindTexture(GL_TEXTURE_2D, TxtIndex);
        }
        glBegin(GL_TRIANGLES);
        for (size_t i=0;i<m.face.size();i++)
        {
            if (!colorPerVert)
                vcg::glColor(m.face[i].C());
            CoordType Pos0=UVTo3DPos(m.face[i].WT(0).P());
            CoordType Pos1=UVTo3DPos(m.face[i].WT(1).P());
            CoordType Pos2=UVTo3DPos(m.face[i].WT(2).P());

            if (TxtIndex>=0)
                vcg::glTexCoord(m.face[i].WT(0).P());
            if (colorPerVert)
                vcg::glColor(m.face[i].V(0)->C());

            vcg::glVertex(Pos0);
            if (TxtIndex>=0)
                vcg::glTexCoord(m.face[i].WT(1).P());
            if (colorPerVert)
                vcg::glColor(m.face[i].V(1)->C());

            vcg::glVertex(Pos1);
            if (TxtIndex>=0)
                vcg::glTexCoord(m.face[i].WT(2).P());
            if (colorPerVert)
                vcg::glColor(m.face[i].V(2)->C());

            vcg::glVertex(Pos2);
        }
        glEnd();
        glPopMatrix();
        glPopAttrib();
    }

};

#endif
