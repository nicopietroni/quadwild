#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <stdio.h>

using namespace vcg;
using namespace std;

class MyPolyFace;
class MyPolyVertex;

struct PolyUsedTypes: public vcg::UsedTypes<vcg::Use<MyPolyVertex>::AsVertexType,
                                            vcg::Use<MyPolyFace>::AsFaceType>{};

class MyPolyVertex:public vcg::Vertex<	PolyUsedTypes,
                                        vcg::vertex::Coord3d,
                                        vcg::vertex::Normal3d,
                                        //vcg::vertex::VFAdj,
                                        vcg::vertex::Mark,
                                        vcg::vertex::BitFlags,
                                        vcg::vertex::CurvatureDird,
                                        vcg::vertex::Qualityd,
                                        vcg::vertex::TexCoord2d>
{};

class MyPolyFace:public vcg::Face<PolyUsedTypes,
                                    vcg::face::PolyInfo,
                                    vcg::face::PFVAdj,
                                    vcg::face::PFFAdj,
                                    vcg::face::BitFlags,
                                    vcg::face::Normal3d,
                                    vcg::face::Color4b,
                                    vcg::face::CurvatureDird,
                                    vcg::face::Qualityd,
                                    vcg::face::WedgeTexCoord2d>
{};


class MyPolyMesh: public vcg::tri::TriMesh< std::vector<MyPolyVertex>,
                                            std::vector<MyPolyFace > >
{
    typedef typename vcg::tri::TriMesh< std::vector<MyPolyVertex>,std::vector<MyPolyFace > > SuperMeshType;

    void UpdateNormal()
    {
        vcg::PolygonalAlgorithm<MyPolyMesh>::UpdateFaceNormalByFitting(*this);
        vcg::tri::UpdateNormal<MyPolyMesh>::PerVertexNormalized(*this);
    }

public:


//    void GLDraw()
//    {

//        glPushAttrib(GL_ALL_ATTRIB_BITS);

//        glDepthRange(0.000001,1);

//        glEnable(GL_LIGHTING);

//        glDisable(GL_CULL_FACE);
//        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

//        vcg::glColor(vcg::Color4b(200,200,200,255));

//        for(unsigned int i=0; i<face.size(); i++)
//        {
//            //vcg::glColor(face[i].C());

//            if(face[i].IsD())  continue;

//            glBegin(GL_POLYGON);

//            vcg::glNormal(face[i].N());
//            for(int j=0; j<face[i].VN(); j++)
//                vcg::glVertex( face[i].V(j)->P() );

//            glEnd();
//        }

//        glDepthRange(0.0,0.999999);
//        glLineWidth(1);
//        glDisable(GL_LIGHTING);
//        vcg::glColor(vcg::Color4b(0,0,0,255));
//        for(unsigned int i=0; i<face.size(); i++)
//        {
//            if(face[i].IsD())  continue;
//            int size=face[i].VN();
//            for(int j=0; j<face[i].VN(); j++)
//            {
//                glBegin(GL_LINES);
//                vcg::glVertex( face[i].V(j)->P() );
//                vcg::glVertex( face[i].V((j+1)%size)->P() );
//                glEnd();
//            }
//        }
//        glPopAttrib();
//    }

    void UpdateAttributesPoly()
    {
        UpdateNormal();
        vcg::tri::UpdateBounding<MyPolyMesh>::Box(*this);
        vcg::tri::UpdateTopology<MyPolyMesh>::FaceFace(*this);
    }

};

MyPolyMesh trimesh,quadmesh;

std::vector<std::vector<size_t> > partitions;
std::vector<std::vector<size_t> > corners;

void LoadPatches(std::string PathName)
{
    FILE *f=fopen(PathName.c_str(),"rt");
    if (f==NULL)
    {
        std::cout<<"ERROR LOADING PATCH FILE"<<std::endl;
        exit(0);
    }
    int NumF;
    fscanf(f,"%d\n",&NumF);
    assert(NumF==trimesh.face.size());
    int MaxPartI=0;
    std::vector<int> FacePar(NumF,-1);

    for (size_t i=0;i<NumF;i++)
    {
        fscanf(f,"%d\n",&FacePar[i]);
        MaxPartI=std::max(FacePar[i],MaxPartI);
    }
    fclose(f);

    partitions.clear();
    partitions.resize(MaxPartI+1);
    for(size_t i=0;i<FacePar.size();i++)
    {
        int CurrP=FacePar[i];
        assert(CurrP>=0);
        assert(CurrP<partitions.size());
        partitions[CurrP].push_back(i);
    }
}

void LoadCorners(std::string PathName)
{
    FILE *f=fopen(PathName.c_str(),"rt");
    if (f==NULL)
    {
        std::cout<<"ERROR LOADING CORNER FILE"<<std::endl;
        exit(0);
    }

    int NumPartitions;
    fscanf(f,"%d\n",&NumPartitions);
    assert(NumPartitions==partitions.size());
    corners.clear();
    corners.resize(NumPartitions);
    for (size_t i=0;i<NumPartitions;i++)
    {
        int NumC;
        fscanf(f,"%d\n",&NumC);
        corners[i].resize(NumC);
        for (size_t j=0;j<NumC;j++)
            fscanf(f,"%d\n",&corners[i][j]);
    }
    fclose(f);
}

//A SIMPLE TEST PROGRAM FOR THE GRID STRUCTURE
int main(int argc, char *argv[])
{

    if(argc<2)
    {
        printf("error: pass one mesh as parameter \n");
        fflush(stdout);
        exit(0);
    }

    //MESH LOAD
    std::string pathM=std::string(argv[1]);
    int mask;
    vcg::tri::io::ImporterOBJ<MyPolyMesh>::LoadMask(pathM.c_str(),mask);
    int err=vcg::tri::io::ImporterOBJ<MyPolyMesh>::Open(trimesh,pathM.c_str(),mask);
    if ((err!=0)&&(err!=5))
    {
        std::cout<<"ERROR LOADING MESH"<<std::endl;
        exit(0);
    }
    std::cout<<"Loaded "<<trimesh.vert.size()<<" vertices"<<std::endl;
    std::cout<<"Loaded "<<trimesh.face.size()<<" faces"<<std::endl;

    //FACE PARTITIONS
    std::string pathP=pathM;
    pathP.erase(pathP.find_last_of("."));
    pathP.append(".patch");
    LoadPatches(pathP);
    std::cout<<"Loaded "<<partitions.size()<<" patches"<<std::endl;

    //PATCH CORNERS
    std::string pathC=pathM;
    pathC.erase(pathC.find_last_of("."));
    pathC.append(".corners");
    LoadCorners(pathC);
    std::cout<<"Loaded "<<corners.size()<<" corners set"<<std::endl;

    //TEST COLORING
    vcg::tri::UpdateColor<MyPolyMesh>::PerFaceConstant(trimesh);
    for(size_t i=0;i<partitions.size();i++)
    {
        vcg::Color4b PartCol=vcg::Color4b::Scatter(partitions.size(),i);
        for(size_t j=0;j<partitions[i].size();j++)
        {
            size_t IndexF=partitions[i][j];
            trimesh.face[IndexF].C()=PartCol;
        }
    }
    vcg::tri::io::ExporterPLY<MyPolyMesh>::Save(trimesh,"out_test.ply",vcg::tri::io::Mask::IOM_FACECOLOR);

}
