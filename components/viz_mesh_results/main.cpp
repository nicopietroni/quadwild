#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <stdio.h>

#include "mesh_types.h"
#include <QFileInfo>
#include <QString>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

std::string pathM,pathF,pathS,pathPatch,pathMeshPatch,pathProject;

TriangleMesh mesh,patch_mesh;

int main(int argc, char *argv[])
{
    //MESH LOAD
    pathM=std::string(argv[1]);
    QString pathMQ=QString(pathM.c_str());
    QFileInfo f_infoM(pathMQ);
    if (!f_infoM.exists())
    {
        std::cout<<"error: mesh fileneme wrong"<<std::endl;
        fflush(stdout);
        exit(0);
    }
    else
        std::cout<<"Mesh file correct"<<std::endl;
    mesh.LoadMesh(pathM.c_str());
    mesh.UpdateAttributes();

    std::cout<<"Loaded "<<mesh.fn<<" faces"<<std::endl;
    std::cout<<"Loaded "<<mesh.vn<<" vertices"<<std::endl;

    pathProject=pathM;
    pathProject.erase(pathProject.find_last_of("."));

    //FIELD LOAD
    pathF=pathProject;
    pathF.append(".rosy");
    QString pathFQ=QString(pathF.c_str());
    QFileInfo f_infoF(pathFQ);
    if (!f_infoF.exists())
    {
        printf("error: field fileneme wrong\n");
        fflush(stdout);
        exit(0);
    }
    else
        std::cout<<"Field file correct"<<std::endl;

    mesh.LoadField(pathF.c_str());
    TriangleMesh FieldMesh;
    mesh.GenerateFieldMesh(FieldMesh);
    std::string pathFMesh=pathProject;
    pathFMesh.append("_field_mesh.ply");
    vcg::tri::io::ExporterPLY<TriangleMesh>::Save(FieldMesh,pathFMesh.c_str(),
                                                  vcg::tri::io::Mask::IOM_FACECOLOR);

    TriangleMesh SingMesh;
    mesh.GenerateSingMesh(SingMesh);
    std::string pathSingMesh=pathProject;
    pathSingMesh.append("_sing_mesh.ply");
    vcg::tri::io::ExporterPLY<TriangleMesh>::Save(SingMesh,pathSingMesh.c_str(),
                                                  vcg::tri::io::Mask::IOM_FACECOLOR);
//    GenerateSingMesh(TriangleMesh &SingMesh,
//                              ScalarType scale=0.01)
    //SHARP LOAD
    pathS=pathProject;
    pathS.append(".sharp");
    QString pathSQ=QString(pathS.c_str());
    QFileInfo f_infoS(pathSQ);
    if (!f_infoS.exists())
    {
        printf("no feature line \n");
        fflush(stdout);
        exit(0);
    }
    else
    {
       std::cout<<"Sharp file correct"<<std::endl;
    }
    mesh.LoadSharpFeatures(pathS);
    TriangleMesh SharpMesh;
    mesh.GenerateEdgeSelMesh(SharpMesh,0.02);
    std::string pathSharpMesh=pathProject;
    pathSharpMesh.append("_sharp_mesh.ply");
    vcg::tri::io::ExporterPLY<TriangleMesh>::Save(SharpMesh,pathSharpMesh.c_str(),
                                                  vcg::tri::io::Mask::IOM_FACECOLOR);


    //PATCH MESH LOAD
    pathMeshPatch=pathProject;
    pathMeshPatch.append("_p0.obj");
    QString pathMPQ=QString(pathMeshPatch.c_str());
    QFileInfo f_infoMP(pathMPQ);
    if (!f_infoMP.exists())
    {
        printf("no mesh patch file \n");
        fflush(stdout);
        exit(0);
    }
    else
    {
       std::cout<<"Patch mesh file correct"<<std::endl;
    }
    patch_mesh.LoadMesh(pathMeshPatch.c_str());
    patch_mesh.UpdateAttributes();
    std::cout<<"Loaded "<<patch_mesh.fn<<" faces"<<std::endl;
    std::cout<<"Loaded "<<patch_mesh.vn<<" vertices"<<std::endl;

    //PATCH LOAD
    pathPatch=pathProject;
    pathPatch.append("_p0.patch");
    QString pathPQ=QString(pathPatch.c_str());
    QFileInfo f_infoP(pathPQ);
    if (!f_infoP.exists())
    {
        printf("no patch file \n");
        fflush(stdout);
        exit(0);
    }
    else
    {
       std::cout<<"Patch file correct"<<std::endl;
    }
    patch_mesh.loadPatchesIntoQ(pathPatch.c_str());
    patch_mesh.SelectPatchBorders();
    patch_mesh.ColorByPartition();

    TriangleMesh BorderPatchMesh;
    patch_mesh.GenerateEdgeSelMesh(BorderPatchMesh,0.02);
    std::string pathBorderPatchMesh=pathProject;
    pathBorderPatchMesh.append("_borderpatch_mesh.ply");
    vcg::tri::io::ExporterPLY<TriangleMesh>::Save(BorderPatchMesh,pathBorderPatchMesh.c_str(),
                                                  vcg::tri::io::Mask::IOM_FACECOLOR);

    std::string ColorPatchMesh=pathProject;
    ColorPatchMesh.append("_colorpatch_mesh.ply");
    vcg::tri::io::ExporterPLY<TriangleMesh>::Save(patch_mesh,ColorPatchMesh.c_str(),
                                                  vcg::tri::io::Mask::IOM_FACECOLOR);


}
