#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <stdio.h>

#include "load_save.h"
#include "mesh_types.h"

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include "quad_from_patches.h"

int main(int argc, char *argv[])
{
    TriangleMesh trimesh;
    std::vector<std::vector<size_t>> trimeshPartitions;
    std::vector<std::vector<size_t>> trimeshCorners;

    PolyMesh quadmesh;
    std::vector<std::vector<size_t>> quadmeshPartitions;
    std::vector<std::vector<size_t>> quadmeshCorners;
    std::vector<int> ilpResult;

    qfp::Parameters parameters;
    parameters.alpha = 0.01;
//    parameters.ilpMethod = qfp::ILPMethod::ABS;

    if(argc<2)
    {
        printf("error: pass one mesh as parameter \n");
        fflush(stdout);
        exit(0);
    }

    //MESH LOAD
    std::string meshFilename = std::string(argv[1]);
    int mask;
    vcg::tri::io::ImporterOBJ<TriangleMesh>::LoadMask(meshFilename.c_str(), mask);
    int err = vcg::tri::io::ImporterOBJ<TriangleMesh>::Open(trimesh, meshFilename.c_str(), mask);

    if ((err!=0)&&(err!=5))
    {
        std::cout<<"ERROR LOADING MESH"<<std::endl;
        exit(0);
    }

    std::cout<<"Loaded "<<trimesh.vert.size()<<" vertices"<<std::endl;
    std::cout<<"Loaded "<<trimesh.face.size()<<" faces"<<std::endl;

    //FACE PARTITIONS
    std::string partitionFilename = meshFilename;
    partitionFilename.erase(partitionFilename.find_last_of("."));
    partitionFilename.append(".patch");
    trimeshPartitions = loadPatches(partitionFilename);
    std::cout<<"Loaded "<<trimeshPartitions.size()<<" patches"<<std::endl;

    //PATCH CORNERS
    std::string cornerFilename = meshFilename;
    cornerFilename.erase(cornerFilename.find_last_of("."));
    cornerFilename.append(".corners");
    trimeshCorners = loadCorners(cornerFilename);
    std::cout<<"Loaded "<<trimeshCorners.size()<<" corners set"<<std::endl;

    //COMPUTE QUADRANGULATION
    qfp::updateAllMeshAttributes(trimesh);
    const std::vector<double> edgeFactor(trimeshPartitions.size(), trimesh.bbox.Diag() / 100);
    qfp::quadrangulationFromPatches(trimesh, trimeshPartitions, trimeshCorners, edgeFactor, parameters, quadmesh, quadmeshPartitions, quadmeshCorners, ilpResult);

    //COLOR AND SAVE QUADRANGULATION
    vcg::tri::UpdateColor<PolyMesh>::PerFaceConstant(quadmesh);
    for(size_t i = 0; i < quadmeshPartitions.size(); i++)
    {
        vcg::Color4b partitionColor = vcg::Color4b::Scatter(quadmeshPartitions.size(),i);
        for(size_t j = 0; j < quadmeshPartitions[i].size(); j++)
        {
            size_t fId = quadmeshPartitions[i][j];
            quadmesh.face[fId].C() = partitionColor;
        }
    }
    vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, "quadrangulation.obj", vcg::tri::io::Mask::IOM_FACECOLOR);


}
