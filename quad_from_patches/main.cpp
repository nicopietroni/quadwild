#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <stdio.h>

#include "load_save.h"
#include "mesh_types.h"

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include "smooth_mesh.h"
#include "quad_from_patches.h"


typename TriangleMesh::ScalarType avgEdge(const TriangleMesh& trimesh);

void loadSetupFile(const std::string& path, qfp::Parameters& parameters, float& scaleFactor);

int main(int argc, char *argv[])
{
    TriangleMesh trimesh;
    std::vector<std::vector<size_t>> trimeshPartitions;
    std::vector<std::vector<size_t>> trimeshCorners;
    std::vector<std::pair<size_t,size_t> > trimeshFeatures;

    PolyMesh quadmesh;
    std::vector<std::vector<size_t>> quadmeshPartitions;
    std::vector<std::vector<size_t>> quadmeshCorners;
    std::vector<int> ilpResult;

//    qfp::Parameters parameters;
//    parameters.alpha = 0.05; //Alpha: blends between isometry (alpha) and regularity (1-alpha)
//    parameters.ilpMethod = qfp::ILPMethod::LEASTSQUARES; //ILP method
//    parameters.timeLimit = 5 * 60; //Time limit in seconds
//    parameters.gapLimit = 0.1; //When it reaches this gap value, optimization stops
//    parameters.minimumGap = 0.25; //Optimization has to reach at least this minimum gap, otherwise faster methods are performed
//    parameters.isometry = true; //Activate isometry
//    parameters.regularityForQuadrilaterals = true; //Activate regularity for quadrilaterals
//    parameters.regularityForNonQuadrilaterals = true; //Activate regularity for non-quadrilaterals
//    parameters.nonQuadrilateralSimilarityFactor = 1.2; //Similarity factor to match sides on non-quad patches
//    parameters.hardParityConstraint = false; //Flag to choose if use hard constraints or not

    qfp::Parameters parameters;
    float scaleFactor;
    loadSetupFile(std::string("basic_setup.txt"), parameters, scaleFactor);

    parameters.hardParityConstraint=true;
    parameters.chartSmoothingIterations = 0;
    parameters.quadrangulationSmoothingIterations = 0; //Fixed borders of the patches

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

    //FEATURES
    std::string featureFilename = meshFilename;
    featureFilename.erase(featureFilename.find_last_of("."));
    featureFilename.append(".feature");
    trimeshFeatures = LoadFeatures(featureFilename);
    std::cout<<"Loaded "<<trimeshFeatures.size()<<" features"<<std::endl;

    //COMPUTE QUADRANGULATION
    qfp::updateAllMeshAttributes(trimesh);
    double EdgeSize=avgEdge(trimesh)*scaleFactor;
    const std::vector<double> edgeFactor(trimeshPartitions.size(), EdgeSize);
    qfp::quadrangulationFromPatches(trimesh, trimeshPartitions, trimeshCorners, edgeFactor, parameters, quadmesh, quadmeshPartitions, quadmeshCorners, ilpResult);

    //COLOR AND SAVE QUADRANGULATION
    vcg::tri::UpdateColor<PolyMesh>::PerFaceConstant(quadmesh);
    for(size_t i = 0; i < quadmeshPartitions.size(); i++)
    {
        vcg::Color4b partitionColor = vcg::Color4b::Scatter(static_cast<int>(quadmeshPartitions.size()), static_cast<int>(i));
        for(size_t j = 0; j < quadmeshPartitions[i].size(); j++)
        {
            size_t fId = quadmeshPartitions[i][j];
            quadmesh.face[fId].C() = partitionColor;
        }
    }

    //SAVE OUTPUT
    std::string outputFilename = meshFilename;
    outputFilename.erase(partitionFilename.find_last_of("."));
    outputFilename.append("_quadrangulation.obj");
    vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, outputFilename.c_str(), vcg::tri::io::Mask::IOM_FACECOLOR);

    //SMOOTH
    std::vector<size_t> QuadPart(quadmesh.face.size(),0);
    for (size_t i=0;i<quadmeshPartitions.size();i++)
        for (size_t j=0;j<quadmeshPartitions[i].size();j++)
            QuadPart[quadmeshPartitions[i][j]]=i;

    std::vector<size_t> TriPart(trimesh.face.size(),0);
    for (size_t i=0;i<trimeshPartitions.size();i++)
        for (size_t j=0;j<trimeshPartitions[i].size();j++)
            TriPart[trimeshPartitions[i][j]]=i;

    SmoothWithFeatures(trimesh,quadmesh,trimeshFeatures,TriPart,QuadPart,Laplacian,10,0.5,EdgeSize);

    //SmoothWithFeatures(trimesh,quadmesh,trimeshFeatures,TriPart,QuadPart,TemplateFit,10,0.5,EdgeSize);

    //SAVE OUTPUT
    outputFilename = meshFilename;
    outputFilename.erase(partitionFilename.find_last_of("."));
    outputFilename.append("_quadrangulation_smooth.obj");
    vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, outputFilename.c_str(), vcg::tri::io::Mask::IOM_FACECOLOR);
}


typename TriangleMesh::ScalarType avgEdge(const TriangleMesh& trimesh)
{
    typedef typename TriangleMesh::ScalarType ScalarType;
    ScalarType AvgVal=0;
    size_t Num=0;
    for (size_t i=0;i<trimesh.face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            AvgVal+=(trimesh.face[i].cP0(j)-trimesh.face[i].cP1(j)).Norm();
            Num++;
        }
    return (AvgVal/Num);
}

void loadSetupFile(const std::string& path, qfp::Parameters& parameters, float& scaleFactor)
{
    FILE *f=fopen(path.c_str(),"rt");
    assert(f!=NULL);

    float alphaF;
    fscanf(f,"alpha %f\n",&alphaF);
    parameters.alpha=alphaF;

    int IntVar=0;
    fscanf(f,"ilpMethod %d\n",&IntVar);
    if (IntVar==0)
        parameters.ilpMethod=qfp::ILPMethod::ABS;
    else
        parameters.ilpMethod=qfp::ILPMethod::LEASTSQUARES;

    float limitF;
    fscanf(f,"timeLimit %f\n",&limitF);
    parameters.timeLimit=limitF;

    float gapF;
    fscanf(f,"gapLimit %f\n",&gapF);
    parameters.gapLimit=gapF;

    float mingapF;
    fscanf(f,"minimumGap %f\n",&mingapF);
    parameters.minimumGap=mingapF;

    IntVar=0;
    fscanf(f,"isometry %d\n",&IntVar);
    if (IntVar==0)
        parameters.isometry=false;
    else
        parameters.isometry=true;

    IntVar=0;
    fscanf(f,"regularityForQuadrilaterals %d\n",&IntVar);
    if (IntVar==0)
        parameters.regularityForQuadrilaterals=false;
    else
        parameters.regularityForQuadrilaterals=true;

    IntVar=0;
    fscanf(f,"regularityForNonQuadrilaterals %d\n",&IntVar);
    if (IntVar==0)
        parameters.regularityForNonQuadrilaterals=false;
    else
        parameters.regularityForNonQuadrilaterals=true;

    float similF;
    fscanf(f,"nonQuadrilateralSimilarityFactor %f\n",&similF);
    parameters.nonQuadrilateralSimilarityFactor=similF;

    IntVar=0;
    fscanf(f,"hardParityConstraint %d\n",&IntVar);
    if (IntVar==0)
        parameters.hardParityConstraint=false;
    else
        parameters.hardParityConstraint=true;

    fscanf(f,"scaleFact %f\n",&scaleFactor);

    fclose(f);
}
