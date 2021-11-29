#ifndef GL_UTILS
#define GL_UTILS


#include <wrap/gl/trimesh.h>

template <class MeshType>
void GLDrawSharpEdges(MeshType &mesh)
{
    typedef typename MeshType::CoordType  CoordType;

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glDepthRange(0,0.9999);
    glLineWidth(5);
    glBegin(GL_LINES);
    for (size_t i=0;i<mesh.face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            if (!mesh.face[i].IsFaceEdgeS(j))continue;

            if (mesh.face[i].FKind[j]==ETConcave)
                vcg::glColor(vcg::Color4b(0,0,255,255));
            else
                vcg::glColor(vcg::Color4b(255,0,0,255));

            CoordType Pos0=mesh.face[i].P0(j);
            CoordType Pos1=mesh.face[i].P1(j);
            vcg::glVertex(Pos0);
            vcg::glVertex(Pos1);
        }
    glEnd();
    glPopAttrib();
}



#endif
