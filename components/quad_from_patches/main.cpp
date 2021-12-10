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

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <stdio.h>
//#include <QFileInfo>
//#include <QDir>

#include "load_save.h"
#include "mesh_types.h"

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include "smooth_mesh.h"
#include "quad_from_patches.h"
#include "quad_mesh_tracer.h"

#include <clocale>

bool LocalUVSm=false;
typename TriangleMesh::ScalarType avgEdge(const TriangleMesh& trimesh);
void loadSetupFile(const std::string& path, QuadRetopology::Parameters& parameters, float& scaleFactor, int& fixedChartClusters);
void SaveSetupFile(const std::string& path, QuadRetopology::Parameters& parameters, float& scaleFactor, int& fixedChartClusters);
//int FindCurrentNum(std::string &pathProject);

int main(int argc, char *argv[])
{
    //Use "." as decimal separator
    std::setlocale(LC_NUMERIC, "en_US.UTF-8");

    int CurrNum=0;

    TriangleMesh trimesh;
    std::vector<std::vector<size_t>> trimeshPartitions;
    std::vector<std::vector<size_t>> trimeshCorners;
    std::vector<std::pair<size_t,size_t> > trimeshFeatures;
    std::vector<size_t> trimeshFeaturesC;

    PolyMesh quadmesh;
    std::vector<std::vector<size_t>> quadmeshPartitions;
    std::vector<std::vector<size_t>> quadmeshCorners;
    std::vector<int> ilpResult;

//    QuadRetopology::Parameters parameters;
//    parameters.alpha = 0.05; //Alpha: blends between isometry (alpha) and regularity (1-alpha)
//    parameters.ilpMethod = QuadRetopology::ILPMethod::LEASTSQUARES; //ILP method
//    parameters.timeLimit = 5 * 60; //Time limit in seconds
//    parameters.gapLimit = 0.1; //When it reaches this gap value, optimization stops
//    parameters.minimumGap = 0.25; //Optimization has to reach at least this minimum gap, otherwise faster methods are performed
//    parameters.isometry = true; //Activate isometry
//    parameters.regularityQuadrilaterals = true; //Activate regularity for quadrilaterals
//    parameters.regularityNonQuadrilaterals = true; //Activate regularity for non-quadrilaterals
//    parameters.regularityNonQuadrilateralWeight = 0.9; //Regularity for non-quadrilaterals weight
//    parameters.alignSingularities = true; //Activate singularity alignment
//    parameters.alignSingularitiesWeight = 0.3; //Singularity alignment weight
//    parameters.repeatLosingConstraintsIterations = true;
//    parameters.repeatLosingConstraintsQuads = false;
//    parameters.repeatLosingConstraintsNonQuads = false;
//    parameters.repeatLosingConstraintsAlign = true;
//    parameters.hardParityConstraint = true; //Flag to choose if use hard constraints or not

    QuadRetopology::Parameters parameters;
    float scaleFactor;
    int fixedChartClusters;
    loadSetupFile(std::string("basic_setup.txt"), parameters, scaleFactor, fixedChartClusters);

    parameters.chartSmoothingIterations = 0; //Chart smoothing
    parameters.quadrangulationFixedSmoothingIterations = 0; //Smoothing with fixed borders of the patches
    parameters.quadrangulationNonFixedSmoothingIterations = 0; //Smoothing with fixed borders of the quadrangulation
    parameters.feasibilityFix = false;

    if(argc<2)
    {
        printf("error: pass one mesh as parameter \n");
        fflush(stdout);
        exit(0);
    }

    if (argc>2)
    {
        CurrNum=atoi(argv[2]);
    }
    //MESH LOAD
    std::string meshFilename = std::string(argv[1]);
    int mask;
    vcg::tri::io::ImporterOBJ<TriangleMesh>::LoadMask(meshFilename.c_str(), mask);
    int err = vcg::tri::io::ImporterOBJ<TriangleMesh>::Open(trimesh, meshFilename.c_str(), mask);
    //trimesh.SolveGeometricIssues();

    if ((err!=0)&&(err!=5))
    {
        std::cout<<"ERROR LOADING MESH"<<std::endl;
        exit(0);
    }
    std::cout<<"MESH NAME "<<meshFilename.c_str()<<std::endl;
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

    //CheckIntegrity(trimesh,trimeshPartitions,trimeshCorners);

    //FEATURES
    std::string featureFilename = meshFilename;
    featureFilename.erase(featureFilename.find_last_of("."));
    featureFilename.append(".feature");
    trimeshFeatures = LoadFeatures(featureFilename);
    std::cout<<"Loaded "<<trimeshFeatures.size()<<" features"<<std::endl;

    //FEATURE CORNERS
    std::string featureCFilename = meshFilename;
    featureCFilename.erase(featureCFilename.find_last_of("."));
    featureCFilename.append(".c_feature");
    trimeshFeaturesC = loadFeatureCorners(featureCFilename);
    std::cout<<"Loaded "<<featureCFilename.size()<<" corner features"<<std::endl;
    loadFeatureCorners(featureCFilename);

    OrientIfNeeded(trimesh,trimeshPartitions,trimeshCorners,trimeshFeatures,trimeshFeaturesC);

    //COMPUTE QUADRANGULATION
    QuadRetopology::internal::updateAllMeshAttributes(trimesh);
    double EdgeSize=avgEdge(trimesh)*scaleFactor;
    std::cout<<"Edge Size "<<EdgeSize<<std::endl;
    const std::vector<double> edgeFactor(trimeshPartitions.size(), EdgeSize);
    qfp::quadrangulationFromPatches(trimesh, trimeshPartitions, trimeshCorners, edgeFactor, parameters, fixedChartClusters, quadmesh, quadmeshPartitions, quadmeshCorners, ilpResult);

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

//    size_t CurrNum=FindCurrentNum(meshFilename);

    //SAVE OUTPUT
    std::string outputFilename = meshFilename;
    outputFilename.erase(partitionFilename.find_last_of("."));
    outputFilename+=std::string("_")+std::to_string(CurrNum)+std::string("_quadrangulation")+std::string(".obj");
    vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, outputFilename.c_str(), vcg::tri::io::Mask::IOM_FACECOLOR);


//    ReMapBoundaries(trimesh,quadmesh,trimeshCorners,trimeshPartitions,
//                    quadmeshCorners,quadmeshPartitions);


    //SMOOTH
    std::vector<size_t> QuadPart(quadmesh.face.size(),0);
    for (size_t i=0;i<quadmeshPartitions.size();i++)
        for (size_t j=0;j<quadmeshPartitions[i].size();j++)
            QuadPart[quadmeshPartitions[i][j]]=i;

    std::vector<size_t> TriPart(trimesh.face.size(),0);
    for (size_t i=0;i<trimeshPartitions.size();i++)
        for (size_t j=0;j<trimeshPartitions[i].size();j++)
            TriPart[trimeshPartitions[i][j]]=i;

    std::vector<size_t> QuadCornersVect;
    for (size_t i=0;i<quadmeshCorners.size();i++)
        for (size_t j=0;j<quadmeshCorners[i].size();j++)
            QuadCornersVect.push_back(quadmeshCorners[i][j]);
    std::sort(QuadCornersVect.begin(),QuadCornersVect.end());
    auto last=std::unique(QuadCornersVect.begin(),QuadCornersVect.end());
    QuadCornersVect.erase(last, QuadCornersVect.end());

    std::cout<<"** SMOOTHING **"<<std::endl;
    //SmoothSubdivide(trimesh,quadmesh,trimeshFeatures,trimeshFeaturesC,TriPart,QuadCornersVect,QuadPart,100,0.5,EdgeSize);
    if (LocalUVSm)
        LocalUVSmooth(quadmesh,trimesh,trimeshFeatures,trimeshFeaturesC,30);
    else
        MultiCostraintSmooth(quadmesh,trimesh,trimeshFeatures,trimeshFeaturesC,TriPart,QuadCornersVect,QuadPart,0.5,EdgeSize,30,1);

    //SAVE OUTPUT
    outputFilename = meshFilename;
    outputFilename.erase(partitionFilename.find_last_of("."));
    //outputFilename.append("_quadrangulation_smooth.obj");
    outputFilename+=std::string("_")+std::to_string(CurrNum)+std::string("_quadrangulation_smooth")+std::string(".obj");

    vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, outputFilename.c_str(), vcg::tri::io::Mask::IOM_FACECOLOR);

#ifdef SAVE_MESHES_FOR_DEBUG
    QuadMeshTracer<PolyMesh> tracerMotorcycle(quadmesh);
    tracerMotorcycle.MotorCycle = true;
    tracerMotorcycle.TracePartitions();
    tracerMotorcycle.SaveColoredMesh("results/quadlayout_motorcycle.obj");
    QuadMeshTracer<PolyMesh> tracer(quadmesh);
    tracer.MotorCycle = false;
    tracer.TracePartitions();
    tracer.SaveColoredMesh("results/quadlayout.obj");
#endif

    for(size_t i=0;i<trimesh.face.size();i++)
        trimesh.face[i].C()=vcg::Color4b::Scatter(trimeshPartitions.size(),TriPart[i]);

#ifdef SAVE_MESHES_FOR_DEBUG
   vcg::tri::io::ExporterOBJ<TriangleMesh>::Save(trimesh,"results/test_tri.obj", vcg::tri::io::Mask::IOM_FACECOLOR);


   std::string setupFilename = meshFilename;
   setupFilename.erase(partitionFilename.find_last_of("."));
   //setupFilename.append("_quadrangulation_setup.txt");
   setupFilename+=std::string("_")+std::to_string(CurrNum)+std::string("_quadrangulation_setup")+std::string(".txt");

   SaveSetupFile(setupFilename, parameters, scaleFactor, fixedChartClusters);
 #endif
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

void loadSetupFile(const std::string& path, QuadRetopology::Parameters& parameters, float& scaleFactor, int& fixedChartClusters)
{
    FILE *f=fopen(path.c_str(),"rt");
    assert(f!=NULL);

    float alphaF;
    fscanf(f,"alpha %f\n",&alphaF);
    std::cout<<"ALPHA "<<alphaF<<std::endl;
    parameters.alpha=alphaF;

    int IntVar=0;
    fscanf(f,"ilpMethod %d\n",&IntVar);
    if (IntVar==0)
        parameters.ilpMethod=QuadRetopology::ILPMethod::ABS;
    else
        parameters.ilpMethod=QuadRetopology::ILPMethod::LEASTSQUARES;

    float limitF;
    fscanf(f,"timeLimit %f\n",&limitF);
    parameters.timeLimit=limitF;

    float gapF;
    fscanf(f,"gapLimit %f\n",&gapF);
    parameters.gapLimit=gapF;

    IntVar=0;
    fscanf(f,"callbackTimeLimit %d",&IntVar);
    parameters.callbackTimeLimit.resize(IntVar);
    for (int i = 0; i < IntVar; i++) {
        fscanf(f," %f", &parameters.callbackTimeLimit[i]);
    }
    fscanf(f,"\n");

    IntVar=0;
    fscanf(f,"callbackGapLimit %d",&IntVar);
    parameters.callbackGapLimit.resize(IntVar);
    for (int i = 0; i < IntVar; i++) {
        fscanf(f," %f", &parameters.callbackGapLimit[i]);
    }
    fscanf(f,"\n");

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
    fscanf(f,"regularityQuadrilaterals %d\n",&IntVar);
    if (IntVar==0)
        parameters.regularityQuadrilaterals=false;
    else
        parameters.regularityQuadrilaterals=true;

    IntVar=0;
    fscanf(f,"regularityNonQuadrilaterals %d\n",&IntVar);
    if (IntVar==0)
        parameters.regularityNonQuadrilaterals=false;
    else
        parameters.regularityNonQuadrilaterals=true;

    float regularityNonQuadrilateralsWeight;
    fscanf(f,"regularityNonQuadrilateralsWeight %f\n",&regularityNonQuadrilateralsWeight);
    parameters.regularityNonQuadrilateralsWeight=regularityNonQuadrilateralsWeight;

    IntVar=0;
    fscanf(f,"alignSingularities %d\n",&IntVar);
    if (IntVar==0)
        parameters.alignSingularities=false;
    else
        parameters.alignSingularities=true;

    float alignSingularitiesWeight;
    fscanf(f,"alignSingularitiesWeight %f\n",&alignSingularitiesWeight);
    parameters.alignSingularitiesWeight=alignSingularitiesWeight;

    IntVar=0;
    fscanf(f,"repeatLosingConstraintsIterations %d\n",&IntVar);
    parameters.repeatLosingConstraintsIterations=IntVar;

    IntVar=0;
    fscanf(f,"repeatLosingConstraintsQuads %d\n",&IntVar);
    if (IntVar==0)
        parameters.repeatLosingConstraintsQuads=false;
    else
        parameters.repeatLosingConstraintsQuads=true;

    IntVar=0;
    fscanf(f,"repeatLosingConstraintsNonQuads %d\n",&IntVar);
    if (IntVar==0)
        parameters.repeatLosingConstraintsNonQuads=false;
    else
        parameters.repeatLosingConstraintsNonQuads=true;

    IntVar=0;
    fscanf(f,"repeatLosingConstraintsAlign %d\n",&IntVar);
    if (IntVar==0)
        parameters.repeatLosingConstraintsAlign=false;
    else
        parameters.repeatLosingConstraintsAlign=true;

    IntVar=0;
    fscanf(f,"hardParityConstraint %d\n",&IntVar);
    if (IntVar==0)
        parameters.hardParityConstraint=false;
    else
        parameters.hardParityConstraint=true;

    fscanf(f,"scaleFact %f\n",&scaleFactor);
    
    fscanf(f,"fixedChartClusters %d\n",&fixedChartClusters);

    fclose(f);
}

void SaveSetupFile(const std::string& path, QuadRetopology::Parameters& parameters, float& scaleFactor, int& fixedChartClusters)
{
    FILE *f=fopen(path.c_str(),"wt");
    assert(f!=NULL);

    fprintf(f,"alpha %f\n", parameters.alpha);

    if (parameters.ilpMethod==QuadRetopology::ILPMethod::ABS)
        fprintf(f,"ilpMethod 0\n");
    else
        fprintf(f,"ilpMethod 1\n");

    fprintf(f,"timeLimit %f\n", parameters.timeLimit);

    fprintf(f,"gapLimit %f\n", parameters.gapLimit);

    fprintf(f,"callbackTimeLimit %d", static_cast<int>(parameters.callbackTimeLimit.size()));
    for (float& time : parameters.callbackTimeLimit) {
        fprintf(f," %f", time);
    }
    fprintf(f,"\n");

    fprintf(f,"callbackGapLimit %d", static_cast<int>(parameters.callbackGapLimit.size()));
    for (float& gap : parameters.callbackGapLimit) {
        fprintf(f," %f", gap);
    }
    fprintf(f,"\n");

    fprintf(f,"minimumGap %f\n", parameters.minimumGap);

    if (parameters.isometry)
        fprintf(f,"isometry 1\n");
    else
        fprintf(f,"isometry 0\n");

    if (parameters.regularityQuadrilaterals)
        fprintf(f,"regularityQuadrilaterals 1\n");
    else
        fprintf(f,"regularityQuadrilaterals 0\n");

    if (parameters.regularityNonQuadrilaterals)
        fprintf(f,"regularityNonQuadrilaterals 1\n");
    else
        fprintf(f,"regularityNonQuadrilaterals 0\n");

    fprintf(f,"regularityNonQuadrilateralsWeight %f\n", parameters.regularityNonQuadrilateralsWeight);

    if (parameters.alignSingularities)
        fprintf(f,"alignSingularities 1\n");
    else
        fprintf(f,"alignSingularities 0\n");

    fprintf(f,"alignSingularitiesWeight %f\n", parameters.alignSingularitiesWeight);

    fprintf(f,"repeatLosingConstraintsIterations %d\n", parameters.repeatLosingConstraintsIterations);

    if (parameters.repeatLosingConstraintsQuads)
        fprintf(f,"repeatLosingConstraintsQuads 1\n");
    else
        fprintf(f,"repeatLosingConstraintsQuads 0\n");

    if (parameters.repeatLosingConstraintsNonQuads)
        fprintf(f,"repeatLosingConstraintsNonQuads 1\n");
    else
        fprintf(f,"repeatLosingConstraintsNonQuads 0\n");

    if (parameters.repeatLosingConstraintsAlign)
        fprintf(f,"repeatLosingConstraintsAlign 1\n");
    else
        fprintf(f,"repeatLosingConstraintsAling 0\n");


    if (parameters.hardParityConstraint)
        fprintf(f,"hardParityConstraint 1\n");
    else
        fprintf(f,"hardParityConstraint 0\n");

    fprintf(f,"scaleFact %f\n", static_cast<float>(scaleFactor));

    fprintf(f,"fixedChartClusters %d\n", fixedChartClusters);

    fclose(f);
}

//int FindCurrentNum(std::string &pathProject)
//{
//    int CurrNum=0;
//    std::string BasePath=pathProject;
//    BasePath.resize(BasePath.find_last_of("/")+1);
//    BasePath="."+BasePath;
//    std::cout<<"basepath "<<BasePath.c_str()<<std::endl;
//    QDir directory(BasePath.c_str());

//    QFile f(pathProject.c_str());
//    QFileInfo fileInfo(f.fileName());
//    QString filename(fileInfo.fileName());
//    std::string NameFile=filename.toStdString();
//    std::cout<<"namefile "<<NameFile.c_str()<<std::endl;
//    NameFile.resize(NameFile.find_last_of("."));
//    std::string Filter=NameFile+"_quadrangulation_smooth*.obj";
//    std::cout<<"filter "<<Filter.c_str()<<std::endl;

////    TestPath+="*.obj";
//  // QStringList projectFiles = directory.entryList(QStringList() <<Filter.c_str(),QDir::Files);
//    QStringList projectFiles = directory.entryList(QStringList(),QDir::Files);
//    CurrNum=0;

//    foreach(QString filename, projectFiles) {
//        int Num;
//        std::string TestScan=NameFile+"%d.obj";
//        sscanf (filename.toStdString().c_str(),TestScan.c_str(),&Num);
//        CurrNum=std::max(CurrNum,(Num+1));
//        std::cout<<"test "<<Num<<std::endl;
//    }
//    return CurrNum;
//}
