//***************************************************************************/
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
#include <iomanip>

#include <triangle_mesh_type.h>
#include <mesh_manager.h>
#include <vcg/space/box3.h>
#include <tracing/mesh_type.h>
#include <tracing/tracer_interface.h>

#include <load_save.h>
#include <mesh_types.h>
#include <smooth_mesh.h>
#include <quad_from_patches.h>
#include <quad_mesh_tracer.h>

#include <clocale>

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

bool DoRemesh=true;
float SharpAngle=35;
float Alpha=0.02;
float ScaleFact=1;
bool HasFeature=false;
std::string PathS;
bool HasField=false;
std::string PathField;

bool loadConfigFile(const std::string & filename)
{

    FILE *f=fopen(filename.c_str(),"rt");

    if (f==NULL)return false;

    std::cout<<"READ CONFIG FILE"<<std::endl;

    int IntVar;
    fscanf(f,"do_remesh %d\n",&IntVar);
    if (IntVar==0)
        DoRemesh=false;
    else
        DoRemesh=true;

    fscanf(f,"sharp_feature_thr %f\n",&SharpAngle);

    fscanf(f,"alpha %f\n",&Alpha);

    fscanf(f,"scaleFact %f\n",&ScaleFact);

    fclose(f);

    std::cout << "[fieldComputation] Successful config import" << std::endl;

    return true;

}

int main(int argc, char *argv[])
{
    //Use "." as decimal separator
    std::setlocale(LC_NUMERIC, "en_US.UTF-8");

    std::cout<<"READING INPUT"<<std::endl;

    if (argc==1)
    {
        std::cout<<"Pass a Mesh (OBJ or PLY) as an argument"<<std::endl;
        exit(0);
    }
    FieldTriMesh tri_mesh;
    bool allQuad;
    std::string PathM=std::string(argv[1]);

    loadConfigFile("basic_setup.txt");

    for (size_t i=2;i<argc;i++)
    {
        int position;
        std::string pathTest=std::string(argv[i]);

        position=pathTest.find(".sharp");
        if (position!=-1)
        {
            PathS=pathTest;
            HasFeature=true;
            continue;
        }

        position=pathTest.find(".txt");
        if (position!=-1)
        {
           loadConfigFile(pathTest.c_str());
           continue;
        }

        position=pathTest.find(".rosy");
        if (position!=-1)
        {
           PathField=pathTest;
           HasField=true;
           continue;
        }
    }


    std::cout<<"Loading:"<<PathM.c_str()<<std::endl;

    bool loaded=tri_mesh.LoadTriMesh(PathM,allQuad);
    tri_mesh.UpdateDataStructures();

    if (!loaded)
    {
        std::cout<<"Mesh Filename Wrong"<<std::endl;
        exit(0);
    }
    std::cout<<"Loaded "<<tri_mesh.fn<<" faces and "<<tri_mesh.vn<<" vertices"<<std::endl;

    std::cout<<"1- REMESH AND FIELD"<<std::endl;

    typename MeshPrepocess<FieldTriMesh>::BatchParam BPar;
    BPar.DoRemesh=DoRemesh;
    BPar.feature_erode_dilate=4;
    BPar.remesher_aspect_ratio=0.3;
    BPar.remesher_iterations=15;
    BPar.remesher_termination_delta=10000;
    BPar.SharpFactor=6;
    BPar.sharp_feature_thr=SharpAngle;
    BPar.surf_dist_check=true;
    BPar.UpdateSharp=(!HasFeature);

    typename vcg::tri::FieldSmoother<FieldTriMesh>::SmoothParam FieldParam;
    FieldParam.alpha_curv=0.3;
    FieldParam.curv_thr=0.8;

    if (HasFeature)
    {
        bool loaded=tri_mesh.LoadSharpFeatures(PathS);
        if (!loaded)
        {
            std::cout<<"ERROR: Wrong Sharp Feature File"<<std::endl;
            exit(0);
        }
        std::cout<<"Sharp Feature Length:"<<tri_mesh.SharpLenght()<<std::endl;
    }
    if (!HasField)
        MeshPrepocess<FieldTriMesh>::BatchProcess(tri_mesh,BPar,FieldParam);
    else
    {
        tri_mesh.LoadField(PathField.c_str());
    }

    MeshPrepocess<FieldTriMesh>::SaveAllData(tri_mesh,PathM);

    std::cout<<"2- TRACING"<<std::endl;

    std::string pathProject=PathM;
    pathProject.erase(pathProject.find_last_of("."));

    PathM=pathProject;
    PathM.append("_rem.obj");
    std::cout<<"Loading Remeshed M:"<<PathM.c_str()<<std::endl;

    std::string PathF=pathProject;
    PathF.append("_rem.rosy");
    std::cout<<"Loading Rosy Field:"<<PathF.c_str()<<std::endl;

    std::string PathS=pathProject;
    PathS.append("_rem.sharp");
    std::cout<<"Loading Sharp F:"<<PathS.c_str()<<std::endl;

    //MESH LOAD
    TraceMesh trace_mesh;
    printf("Loading the mesh \n");
    bool loadedMesh=trace_mesh.LoadMesh(PathM);
    assert(loadedMesh);
    trace_mesh.UpdateAttributes();

    //FIELD LOAD
    bool loadedField=trace_mesh.LoadField(PathF);
    assert(loadedField);
    trace_mesh.UpdateAttributes();

    //SHARP LOAD
    bool loadedFeatures=trace_mesh.LoadSharpFeatures(PathS);
    assert(loadedFeatures);
    trace_mesh.SolveGeometricIssues();
    trace_mesh.UpdateSharpFeaturesFromSelection();

    //preprocessing mesh
    PreProcessMesh(trace_mesh);

    //initializing graph
    VertexFieldGraph<TraceMesh> VGraph(trace_mesh);
    VGraph.InitGraph(false);

    //INIT TRACER
    typedef PatchTracer<TraceMesh> TracerType;
    TracerType PTr(VGraph);
    TraceMesh::ScalarType Drift=100;
    bool add_only_needed=true;
    bool final_removal=true;
    bool meta_mesh_collapse=true;
    bool force_split=false;
    PTr.sample_ratio=0.01;
    PTr.CClarkability=1;
    PTr.split_on_removal=true;
    PTr.away_from_singular=true;
    PTr.match_valence=true;
    PTr.check_quality_functor=false;
    PTr.MinVal=3;
    PTr.MaxVal=5;
    PTr.Concave_Need=1;

    //TRACING
    PTr.InitTracer(Drift,false);
    RecursiveProcess<TracerType>(PTr,Drift, add_only_needed,final_removal,true,meta_mesh_collapse,force_split,true,false);
    PTr.SmoothPatches();
    SaveAllData(PTr,pathProject,0,false,false);


    std::cout<<"3- QUADDING"<<std::endl;
    TriangleMesh to_quad_trimesh;
    std::vector<std::vector<size_t>> trimeshPartitions;
    std::vector<std::vector<size_t>> trimeshCorners;
    std::vector<std::pair<size_t,size_t> > trimeshFeatures;
    std::vector<size_t> trimeshFeaturesC;

    PolyMesh quadmesh;
    std::vector<std::vector<size_t>> quadmeshPartitions;
    std::vector<std::vector<size_t>> quadmeshCorners;
    std::vector<int> ilpResult;

    PathM=pathProject;
    PathM.append("_p0.obj");
    int mask;
    vcg::tri::io::ImporterOBJ<TriangleMesh>::LoadMask(PathM.c_str(), mask);
    int err = vcg::tri::io::ImporterOBJ<TriangleMesh>::Open(to_quad_trimesh, PathM.c_str(), mask);
    if ((err!=0)&&(err!=5))
        assert(0);

    //FACE PARTITIONS
    std::string partitionFilename = pathProject;
    partitionFilename.append("_p0.patch");
    trimeshPartitions = loadPatches(partitionFilename);
    std::cout<<"Loaded "<<trimeshPartitions.size()<<" patches"<<std::endl;

    //PATCH CORNERS
    std::string cornerFilename = pathProject;
    cornerFilename.append("_p0.corners");
    trimeshCorners = loadCorners(cornerFilename);
    std::cout<<"Loaded "<<trimeshCorners.size()<<" corners set"<<std::endl;

    //FEATURES
    std::string featureFilename = pathProject;
    featureFilename.append("_p0.feature");
    trimeshFeatures = LoadFeatures(featureFilename);
    std::cout<<"Loaded "<<trimeshFeatures.size()<<" features"<<std::endl;

    //FEATURE CORNERS
    std::string featureCFilename = pathProject;
    featureCFilename.append("_p0.c_feature");
    trimeshFeaturesC = loadFeatureCorners(featureCFilename);
    std::cout<<"Loaded "<<featureCFilename.size()<<" corner features"<<std::endl;
    loadFeatureCorners(featureCFilename);

    OrientIfNeeded(to_quad_trimesh,trimeshPartitions,trimeshCorners,trimeshFeatures,trimeshFeaturesC);

    //COMPUTE QUADRANGULATION
    QuadRetopology::internal::updateAllMeshAttributes(to_quad_trimesh);

    QuadRetopology::Parameters parameters;
    float scaleFactor;
    int fixedChartClusters;

    parameters.alpha=Alpha;
    parameters.ilpMethod=QuadRetopology::ILPMethod::LEASTSQUARES;
    parameters.timeLimit=200;
    parameters.gapLimit=0.0;
    parameters.callbackTimeLimit.push_back(3.0);
    parameters.callbackTimeLimit.push_back(5.0);
    parameters.callbackTimeLimit.push_back(10.0);
    parameters.callbackTimeLimit.push_back(20.0);
    parameters.callbackTimeLimit.push_back(30.0);
    parameters.callbackTimeLimit.push_back(60.0);
    parameters.callbackTimeLimit.push_back(90.0);
    parameters.callbackTimeLimit.push_back(120.0);

    parameters.callbackGapLimit.push_back(0.005);
    parameters.callbackGapLimit.push_back(0.02);
    parameters.callbackGapLimit.push_back(0.05);
    parameters.callbackGapLimit.push_back(0.1);
    parameters.callbackGapLimit.push_back(0.15);
    parameters.callbackGapLimit.push_back(0.20);
    parameters.callbackGapLimit.push_back(0.25);
    parameters.callbackGapLimit.push_back(0.3);

    parameters.minimumGap=0.4;

    parameters.isometry=true;

    parameters.regularityQuadrilaterals=true;

    parameters.regularityNonQuadrilaterals=true;

    parameters.regularityNonQuadrilateralsWeight=0.9;

    parameters.alignSingularities=true;

    parameters.alignSingularitiesWeight=0.1;

    parameters.repeatLosingConstraintsIterations=true;

    parameters.repeatLosingConstraintsQuads=false;

    parameters.repeatLosingConstraintsNonQuads=false;

    parameters.repeatLosingConstraintsAlign=true;

    parameters.hardParityConstraint=true;

    scaleFactor=ScaleFact;

    fixedChartClusters=300;

    parameters.chartSmoothingIterations = 0; //Chart smoothing
    parameters.quadrangulationFixedSmoothingIterations = 0; //Smoothing with fixed borders of the patches
    parameters.quadrangulationNonFixedSmoothingIterations = 0; //Smoothing with fixed borders of the quadrangulation
    parameters.feasibilityFix = false;

    double EdgeSize=avgEdge(to_quad_trimesh)*scaleFactor;
    std::cout<<"Edge Size "<<EdgeSize<<std::endl;
    const std::vector<double> edgeFactor(trimeshPartitions.size(), EdgeSize);

    qfp::quadrangulationFromPatches(to_quad_trimesh, trimeshPartitions, trimeshCorners, edgeFactor, parameters, fixedChartClusters, quadmesh, quadmeshPartitions, quadmeshCorners, ilpResult);

    //SAVE OUTPUT
    std::string outputFilename = pathProject;
    outputFilename+=std::string("_quadrangulation")+std::string(".obj");
    vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, outputFilename.c_str(),0);


    //SMOOTH
    std::vector<size_t> QuadPart(quadmesh.face.size(),0);
    for (size_t i=0;i<quadmeshPartitions.size();i++)
        for (size_t j=0;j<quadmeshPartitions[i].size();j++)
            QuadPart[quadmeshPartitions[i][j]]=i;

    std::vector<size_t> TriPart(to_quad_trimesh.face.size(),0);
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
    MultiCostraintSmooth(quadmesh,to_quad_trimesh,trimeshFeatures,trimeshFeaturesC,TriPart,QuadCornersVect,QuadPart,0.5,EdgeSize,30,1);

    //SAVE OUTPUT
    outputFilename = pathProject;
    outputFilename+=std::string("_quadrangulation_smooth")+std::string(".obj");

    vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, outputFilename.c_str(),0);

}
