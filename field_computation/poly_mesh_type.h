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

#ifndef POLY_MESH_TYPE
#define POLY_MESH_TYPE

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>

#ifdef MIQ_QUADRANGULATE
class PolyFace;
class PolyVertex;

struct PUsedTypes: public vcg::UsedTypes<vcg::Use<PolyVertex>  ::AsVertexType,
        vcg::Use<PolyFace>	::AsFaceType>{};

class PolyVertex:public vcg::Vertex<	PUsedTypes,
        vcg::vertex::Coord3d,
        vcg::vertex::Normal3d,
        vcg::vertex::Mark,
        vcg::vertex::BitFlags,
        vcg::vertex::Qualityd,
        vcg::vertex::TexCoord2d>{} ;

class PolyFace:public vcg::Face<
        PUsedTypes
        ,vcg::face::PolyInfo
        ,vcg::face::PFVAdj
        ,vcg::face::PFFAdj
        ,vcg::face::BitFlags
        ,vcg::face::Normal3d
        ,vcg::face::Color4b
        ,vcg::face::Qualityd      // face quality.
        ,vcg::face::BitFlags
        ,vcg::face::Mark
        ,vcg::face::CurvatureDird> {
};

class PMesh: public
        vcg::tri::TriMesh<
        std::vector<PolyVertex>,	// the vector of vertices
        std::vector<PolyFace >     // the vector of faces
        >
{
public:

    void TriangulateQuadBySplit(size_t IndexF)
    {

        size_t sizeV=face[IndexF].VN();
        assert(sizeV==4);

        //then reupdate the faces
        VertexType * v0=face[IndexF].V(0);
        VertexType * v1=face[IndexF].V(1);
        VertexType * v2=face[IndexF].V(2);
        VertexType * v3=face[IndexF].V(3);

        face[IndexF].Dealloc();
        face[IndexF].Alloc(3);
        face[IndexF].V(0)=v0;
        face[IndexF].V(1)=v1;
        face[IndexF].V(2)=v3;

        vcg::tri::Allocator<PMesh>::AddFaces(*this,1);

        face.back().Alloc(3);
        face.back().V(0)=v1;
        face.back().V(1)=v2;
        face.back().V(2)=v3;
    }

    void TriangulateQuadBySplit()
    {
        for (size_t i=0;i<face.size();i++)
        {
            if(face[i].VN()!=4)continue;
            TriangulateQuadBySplit(i);
        }
    }
    void GLDraw(bool DrawEdges=true)
    {

        glPushAttrib(GL_ALL_ATTRIB_BITS);

        glDepthRange(0.000001,1);

        glEnable(GL_LIGHTING);

        glDisable(GL_CULL_FACE);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

        for(unsigned int i=0; i<face.size(); i++)
        {
            vcg::glColor(face[i].C());

            if(face[i].IsD())  continue;
            //if(face[i].Filtered) continue;

            vcg::glNormal(face[i].N());

            glBegin(GL_POLYGON);
            for(int j=0; j<face[i].VN(); j++)
                vcg::glVertex( face[i].V(j)->P() );

            glEnd();

        }

        //DrawEdgeSmoothPairs();

        if (DrawEdges)
        {
            glDepthRange(0.0,0.999999);
            glDisable(GL_LIGHTING);

            //            vcg::glColor(vcg::Color4b(255,0,0,255));
            //            for(unsigned int i=0; i<vert.size(); i++)
            //            {
            //                if (!vert[i].Irregular)continue;
            //                glPointSize(20);
            //                glBegin(GL_POINTS);
            //                vcg::glVertex( vert[i].P());
            //                glEnd();
            //            }
            for(unsigned int i=0; i<face.size(); i++)
            {
                if(face[i].IsD())  continue;
                int size=face[i].VN();
                for(int j=0; j<face[i].VN(); j++)
                {

                    glLineWidth(4);
                    vcg::glColor(vcg::Color4b(0,0,0,255));

                    CoordType pos0=face[i].V(j)->P();
                    CoordType pos1=face[i].V((j+1)%size)->P();

                    glBegin(GL_LINES);
                    vcg::glVertex( pos0);
                    vcg::glVertex( pos1);
                    glEnd();
                }
            }
        }
        //glEnd();
        glPopAttrib();

    }

    void UpdateNormal()
    {
        vcg::PolygonalAlgorithm<PMesh>::UpdateFaceNormalByFitting(*this);
        vcg::tri::UpdateNormal<PMesh>::PerVertexNormalized(*this);
    }

    void UpdateAttributes()
    {
        UpdateNormal();
        vcg::tri::UpdateBounding<PMesh>::Box(*this);
        vcg::tri::UpdateTopology<PMesh>::FaceFace(*this);
        //vcg::PolygonalAlgorithm<QuadMeshC>::UpdateBorderVertexFromPFFAdj(*this);
        vcg::tri::UpdateFlags<PMesh>::VertexBorderFromFaceAdj(*this);

    }
};
#endif



#endif
