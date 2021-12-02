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

#ifndef GL_META_MESH
#define GL_META_MESH

#include "metamesh.h"

template <class MeshType>
class GLMetaMesh
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef MetaMesh<MeshType> MetaMeshType;
    //typedef typename vcg::face::Pos<FaceType> PosType;

public:

    enum CEdgeMode{CMNone,CMRemovable,CMLenght,CMPathId};

    static void GlDrawMetaVert(MetaMeshType &m,size_t IndexV)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9995);

        if (m.MVert[IndexV].FixedEnd)
        {
            glPointSize(20);
            glColor(vcg::Color4b(255,0,0,255));
        }
        else
        {
            glPointSize(10);
            glColor(vcg::Color4b(0,255,0,255));
        }

        CoordType Pos = m.MetaVertPos(IndexV);
        glBegin(GL_POINTS);
        vcg::glVertex(Pos);
        glEnd();
        glPopAttrib();
    }

    static void GlDrawMetaFaceAdj(MetaMeshType &m,size_t IndexF)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9995);

        glLineWidth(5);


        size_t sizeV=m.MFaces[IndexF].V.size();

        CoordType Bary0=m.MFaces[IndexF].BaryF;
        for (size_t i=0;i<sizeV;i++)
        {
            CoordType PosE0= m.MetaVertPos(m.MFaces[IndexF].V[i]);
            CoordType PosE1= m.MetaVertPos(m.MFaces[IndexF].V[(i+1)%sizeV]);
            CoordType AvgE=(PosE0+PosE1)/2;
            int AdjF=m.MFaces[IndexF].AdjF[i].first;

            CoordType Pos0=Bary0*0.5+AvgE*0.5;
            CoordType Pos1=AvgE;
            CoordType Pos2=AvgE;
            vcg::Color4b adj_col(255,0,255,255);
            if (AdjF!=-1)
            {
                assert(m.MFaces[IndexF].AdjF[i].second!=-1);
                int AdjF=m.MFaces[IndexF].AdjF[i].first;
                //int AdjE=MFaces[IndexF].AdjF[i].second;
                CoordType Bary1=m.MFaces[AdjF].BaryF;
                Pos2=Bary1*0.5+AvgE*0.5;//Bary1;
                adj_col=vcg::Color4b(0,0,255,255);
            }
            glColor(adj_col);
            glBegin(GL_LINE_STRIP);
            vcg::glVertex(Pos0);
            vcg::glVertex(Pos1);
            vcg::glVertex(Pos2);
            glEnd();
        }
        glPopAttrib();
    }

    static void GlDrawMetaFaceVert(MetaMeshType &m,size_t IndexF)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9995);

        glPointSize(10);
        glColor(vcg::Color4b(0,0,0,255));

        size_t sizeV=m.MFaces[IndexF].V.size();

        CoordType bary=m.MFaces[IndexF].BaryF;

        glBegin(GL_POINTS);
        for (size_t i=0;i<sizeV;i++)
        {
            if (!m.MFaces[IndexF].RealV[i])continue;
            CoordType pos=m.MetaVertPos(m.MFaces[IndexF].V[i]);
            pos=pos*0.8+bary*0.2;
            vcg::glVertex(pos);
        }
        glEnd();


        glPopAttrib();
    }

    static void GlDrawMetaFaceVertField(MetaMeshType &m,size_t IndexF)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9995);


        size_t sizeV=m.MFaces[IndexF].V.size();

        CoordType bary=m.MFaces[IndexF].BaryF;

        for (size_t i=0;i<sizeV;i++)
        {
            CoordType pos0=m.MetaVertPos(m.MFaces[IndexF].V[i]);
            pos0=pos0*0.8+bary*0.2;

            std::vector<CoordType> EdgeDir;
            m.GetEdgeDirOnV(IndexF,i,EdgeDir);
            assert(EdgeDir.size()==2);
            ScalarType sizeD=(*m.mesh).bbox.Diag()*0.01;
            CoordType pos1=pos0+EdgeDir[0]*sizeD;
            CoordType pos2=pos0+EdgeDir[1]*sizeD;
            if (m.IsFieldCornerFaceV(IndexF,i))
            {
               glLineWidth(10);
               glColor(vcg::Color4b(255,255,0,255));
            }
            else
            {
               glLineWidth(5);
               glColor(vcg::Color4b(200,200,200,255));
            }
            glBegin(GL_LINES);
            vcg::glVertex(pos0);
            vcg::glVertex(pos1);
            vcg::glVertex(pos0);
            vcg::glVertex(pos2);
            glEnd();
        }



        glPopAttrib();
    }

    static void GlDrawMetaFaceSides(MetaMeshType &m,size_t IndexF)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9995);
        glDisable(GL_CULL_FACE);
        //glLineWidth(5);


        size_t sizeV=m.NumSides(IndexF);
        int ExpVal=m.MFaces[IndexF].ExpectedVal;
        if ((sizeV<3)||(sizeV>6))
            glColor(vcg::Color4b(255,0,0,255));
        if (sizeV==3)
            glColor(vcg::Color4b(0,255,255,255));
        if (sizeV==4)
            glColor(vcg::Color4b(100,100,100,255));
        if (sizeV==5)
            glColor(vcg::Color4b(255,0,255,255));
        if (sizeV==6)
            glColor(vcg::Color4b(255,255,0,255));
        if (ExpVal!=sizeV)
            glColor(vcg::Color4b(255,0,0,255));

        CoordType bary=m.MFaces[IndexF].BaryF;

        //        for (size_t i=0;i<sizeV;i++)
        //        {
        //            size_t IndexV0,IndexV1;
        //            GetSideExtremes(IndexF,i,IndexV0,IndexV1);

        //            CoordType Pos0= MetaVertPos(MFaces[IndexF].V[IndexV0]);
        //            Pos0=bary*0.75+Pos0*0.25;
        //            CoordType Pos1= MetaVertPos(MFaces[IndexF].V[IndexV1]);
        //            Pos1=bary*0.75+Pos1*0.25;
        //            glBegin(GL_LINES);
        //            vcg::glVertex(Pos0);
        //            vcg::glVertex(Pos1);
        //            glEnd();
        //        }

        glBegin(GL_POLYGON);
        for (size_t i=0;i<sizeV;i++)
        {
            size_t IndexV0,IndexV1;
            m.GetSideExtremes(IndexF,i,IndexV0,IndexV1);

            CoordType Pos0= m.MetaVertPos(m.MFaces[IndexF].V[IndexV0]);
            Pos0=bary*0.75+Pos0*0.25;

            vcg::glVertex(Pos0);
        }
        glEnd();

        glPopAttrib();
    }



    static void GlDrawMetaFace(MetaMeshType &m,size_t IndexF,const CEdgeMode EMode)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9995);

        glLineWidth(5);


        size_t sizeV=m.MFaces[IndexF].V.size();
        for (size_t i=0;i<sizeV;i++)
        {
            if (EMode==CMNone)
                glColor(vcg::Color4b(200,200,200,255));
            if (EMode==CMLenght)
                glColor(vcg::Color4b::ColorRamp(0,m.MaxL,m.MFaces[IndexF].MetaL[i]));
            if (EMode==CMRemovable)
            {
                if (m.IsRemovableEdgeVisual(IndexF,i))
                    glColor(vcg::Color4b(0,255,0,255));
                else
                    glColor(vcg::Color4b(255,0,0,255));
            }
            if (EMode==CMPathId)
            {
                size_t IndexP=m.MFaces[IndexF].PathId[i];
                glColor(vcg::Color4b::Scatter(m.MaxPath,IndexP));
            }
            CoordType Pos0= m.MetaVertPos(m.MFaces[IndexF].V[i]);
            CoordType Pos1= m.MetaVertPos(m.MFaces[IndexF].V[(i+1)%sizeV]);
            glBegin(GL_LINES);
            vcg::glVertex(Pos0);
            vcg::glVertex(Pos1);
            glEnd();
        }
        glPopAttrib();
    }


    static void GLDraw(MetaMeshType &m,
                       const CEdgeMode EMode=CMRemovable,
                       bool DrawMetaEdgeAdj=true)
    {
        for (size_t i=0;i<m.MVert.size();i++)
            GlDrawMetaVert(m,i);

        for (size_t i=0;i<m.MFaces.size();i++)
        {
            GlDrawMetaFace(m,i,EMode);
            GlDrawMetaFaceSides(m,i);

            if (DrawMetaEdgeAdj)
                GlDrawMetaFaceAdj(m,i);

            GlDrawMetaFaceVert(m,i);
            GlDrawMetaFaceVertField(m,i);
        }
        //        for (size_t i=0;i<MFaces.size();i++)
        //            GlDrawMetaFace(i);
    }

};
#endif
