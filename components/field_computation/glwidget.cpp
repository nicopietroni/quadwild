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

#include <GL/glew.h>
#include <QMouseEvent>

#include <chrono>

#include <math.h>
#include "glwidget.h"
#include "mesh_manager.h"
#include "fields/field_smoother.h"
#include "triangle_mesh_type.h"
#include "poly_mesh_type.h"
#include <wrap/qt/trackball.h>
#include <wrap/gl/picking.h>
#include <wrap/qt/anttweakbarMapper.h>
#include <wrap/io_trimesh/import_field.h>
#include <wrap/io_trimesh/export_field.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/gl/trimesh.h>
//#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
#include <wrap/gl/gl_field.h>
#include <AutoRemesher.h>
#include <vcg/complex/algorithms/hole.h>
#include "mesh_field_smoother.h"
#include "gl_utils.h"

#ifdef MIQ_QUADRANGULATE
#include <wrap/igl/miq_parametrization.h>
#include <vcg/complex/algorithms/quadrangulator.h>
#endif

std::string pathM="";
std::string pathS="";

vcg::Trackball track;//the active manipulator

bool drawfield=false;

FieldTriMesh tri_mesh;

#ifdef MIQ_QUADRANGULATE
PMesh quad_mesh;
bool quadrangulated=false;
vcg::tri::MiQParametrizer<FieldTriMesh>::MIQParameters MiqP;
#endif

TwBar *barQuad;

vcg::GlTrimesh<FieldTriMesh> glWrap;
vcg::GLField<FieldTriMesh> glField;

typedef typename FieldTriMesh::ScalarType ScalarType;
typedef typename FieldTriMesh::CoordType CoordType;

int Iterations;
ScalarType EdgeStep;
ScalarType Multiplier=2;
ScalarType SharpFactor=6;
ScalarType alpha=0.3;
//ScalarType sharp_feature_thr=45;


int xMouse,yMouse;

//vcg::GridStaticPtr<FieldTriMesh::FaceType,FieldTriMesh::ScalarType> Gr;

typedef vcg::tri::FieldSmoother<FieldTriMesh> FieldSmootherType;
FieldSmootherType::SmoothParam FieldParam;

AutoRemesher<FieldTriMesh>::Params RemPar;
bool do_batch=false;
bool has_features=false;
bool has_features_fl=false;
//size_t ErodeDilateSteps=4;

bool do_remesh=true;
int remesher_iterations=15;
ScalarType remesher_aspect_ratio=0.3;
int remesher_termination_delta=10000;
bool surf_dist_check=true;
ScalarType sharp_feature_thr=35;
int feature_erode_dilate=4;

//MeshPrepocess<FieldTriMesh> MP(tri_mesh);

void InitSharp()
{
    assert(!(has_features || has_features_fl));
    tri_mesh.UpdateDataStructures();
    tri_mesh.InitSharpFeatures(sharp_feature_thr);
    //MP.InitSharpFeatures(sharp_feature_thr);
}


//void DoAutoRemesh()
//{
//    RemPar.iterations   = remesher_iterations;
//    RemPar.targetAspect = remesher_aspect_ratio;
//    RemPar.targetDeltaFN= remesher_termination_delta;
//    RemPar.creaseAngle  = sharp_feature_thr;
//    RemPar.erodeDilate  = feature_erode_dilate;
//    RemPar.userSelectedCreases = true;
//    RemPar.surfDistCheck = true;

//    auto t1 = std::chrono::high_resolution_clock::now();
//    std::cout << "cleaning the mesh..." << std::endl;
//    std::shared_ptr<FieldTriMesh> clean = AutoRemesher<FieldTriMesh>::CleanMesh(tri_mesh, false);
//    auto t2 = std::chrono::high_resolution_clock::now();
//    std::cout << "Cleaning time: " << std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << std::endl;

//    t1 = std::chrono::high_resolution_clock::now();
//    clean->UpdateDataStructures();

////    MP.InitSharpFeatures(sharp_feature_thr);
////    MP.ErodeDilate(feature_erode_dilate);
//    clean->InitSharpFeatures(sharp_feature_thr);
//    clean->ErodeDilate(feature_erode_dilate);

//    t2 = std::chrono::high_resolution_clock::now();
//    std::cout << "Sharp Features time: " << std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << std::endl;

//    t1 = std::chrono::high_resolution_clock::now();
//    std::shared_ptr<FieldTriMesh> ret=AutoRemesher<FieldTriMesh>::Remesh(*clean,RemPar);
//    t2 = std::chrono::high_resolution_clock::now();
//    std::cout << "Remesh time: " << std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << std::endl;
//    AutoRemesher<FieldTriMesh>::SplitNonManifold(*ret);

//    tri_mesh.Clear();
//    vcg::tri::Append<FieldTriMesh,FieldTriMesh>::Mesh(tri_mesh,(*ret));
//    vcg::tri::Clean<FieldTriMesh>::RemoveUnreferencedVertex(tri_mesh);
//    vcg::tri::Allocator<FieldTriMesh>::CompactEveryVector(tri_mesh);
//    tri_mesh.UpdateDataStructures();
//    //MeshPrepocess<FieldTriMesh>::AutoRemesh(tri_mesh,RemPar);
//    Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
//}

//void DoAutoRemesh()
//{
//    RemPar.iterations   = remesher_iterations;
//    RemPar.targetAspect = remesher_aspect_ratio;
//    RemPar.targetDeltaFN= remesher_termination_delta;
//    RemPar.creaseAngle  = sharp_feature_thr;
//    RemPar.erodeDilate  = feature_erode_dilate;
//    RemPar.userSelectedCreases = true;
//    RemPar.surfDistCheck = false;

//    auto t1 = std::chrono::high_resolution_clock::now();
//    std::cout << "cleaning the mesh..." << std::endl;
//    std::shared_ptr<FieldTriMesh> clean = AutoRemesher<FieldTriMesh>::CleanMesh(tri_mesh, false);
//    auto t2 = std::chrono::high_resolution_clock::now();
//    std::cout << "Cleaning time: " << std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << std::endl;

//    t1 = std::chrono::high_resolution_clock::now();
//    clean->UpdateDataStructures();

////    clean->InitSharpFeatures(sharp_feature_thr);
////    clean->ErodeDilate(feature_erode_dilate);

//    t2 = std::chrono::high_resolution_clock::now();
//    std::cout << "Sharp Features time: " << std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << std::endl;

//    t1 = std::chrono::high_resolution_clock::now();
//    std::shared_ptr<FieldTriMesh> ret=AutoRemesher<FieldTriMesh>::Remesh(*clean,RemPar);
//    t2 = std::chrono::high_resolution_clock::now();
//    std::cout << "Remesh time: " << std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << std::endl;
//    AutoRemesher<FieldTriMesh>::SplitNonManifold(*ret);

//    tri_mesh.Clear();
//    vcg::tri::Append<FieldTriMesh,FieldTriMesh>::Mesh(tri_mesh,(*ret));
//    vcg::tri::Clean<FieldTriMesh>::RemoveUnreferencedVertex(tri_mesh);
//    vcg::tri::Allocator<FieldTriMesh>::CompactEveryVector(tri_mesh);
//    tri_mesh.UpdateDataStructures();

//    tri_mesh.InitSharpFeatures(sharp_feature_thr);
//    tri_mesh.ErodeDilate(feature_erode_dilate);

//    //MeshPrepocess<FieldTriMesh>::AutoRemesh(tri_mesh,RemPar);
//    Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
//}


//void DoAutoRemesh()
//{
//    RemPar.iterations   = remesher_iterations;
//    RemPar.targetAspect = remesher_aspect_ratio;
//    RemPar.targetDeltaFN= remesher_termination_delta;
//    RemPar.creaseAngle  = sharp_feature_thr;
//    RemPar.erodeDilate  = feature_erode_dilate;
//    RemPar.userSelectedCreases = true;
//    RemPar.surfDistCheck = false;

//    auto t1 = std::chrono::high_resolution_clock::now();
//    std::cout << "cleaning the mesh..." << std::endl;
//    std::shared_ptr<FieldTriMesh> clean = AutoRemesher<FieldTriMesh>::CleanMesh(tri_mesh, false);
//    auto t2 = std::chrono::high_resolution_clock::now();
//    std::cout << "Cleaning time: " << std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << std::endl;

//    t1 = std::chrono::high_resolution_clock::now();
//    clean->UpdateDataStructures();

////    clean->InitSharpFeatures(sharp_feature_thr);
////    clean->ErodeDilate(feature_erode_dilate);

//    t2 = std::chrono::high_resolution_clock::now();
//    std::cout << "Sharp Features time: " << std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << std::endl;

//    t1 = std::chrono::high_resolution_clock::now();
//    std::shared_ptr<FieldTriMesh> ret=AutoRemesher<FieldTriMesh>::Remesh(*clean,RemPar);
//    t2 = std::chrono::high_resolution_clock::now();
//    std::cout << "Remesh time: " << std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << std::endl;
//    AutoRemesher<FieldTriMesh>::SplitNonManifold(*ret);

//    tri_mesh.Clear();
//    vcg::tri::Append<FieldTriMesh,FieldTriMesh>::Mesh(tri_mesh,(*ret));
//    vcg::tri::Clean<FieldTriMesh>::RemoveUnreferencedVertex(tri_mesh);
//    vcg::tri::Allocator<FieldTriMesh>::CompactEveryVector(tri_mesh);
//    tri_mesh.UpdateDataStructures();

//    tri_mesh.InitSharpFeatures(sharp_feature_thr);
//    tri_mesh.ErodeDilate(feature_erode_dilate);

//    //MeshPrepocess<FieldTriMesh>::AutoRemesh(tri_mesh,RemPar);
//    Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
//}

void SaveAllData()
{
    MeshPrepocess<FieldTriMesh>::SaveAllData(tri_mesh,pathM);
//    std::string projM=pathM;
//    size_t indexExt=projM.find_last_of(".");
//    projM=projM.substr(0,indexExt);
//    std::string meshName=projM+std::string("_rem.obj");
//    std::string fieldName=projM+std::string("_rem.rosy");
//    std::string sharpName=projM+std::string("_rem.sharp");
//    std::cout<<"Saving Mesh TO:"<<meshName.c_str()<<std::endl;
//    std::cout<<"Saving Field TO:"<<fieldName.c_str()<<std::endl;
//    std::cout<<"Saving Sharp TO:"<<sharpName.c_str()<<std::endl;
//    tri_mesh.SaveTriMesh(meshName.c_str());
//    tri_mesh.SaveField(fieldName.c_str());
//    tri_mesh.SaveSharpFeatures(sharpName.c_str());


////    std::string orFaceName=projM+std::string("_rem_origf.txt");
////    std::cout<<"Saving Field TO:"<<orFaceName.c_str()<<std::endl;
////    tri_mesh.SaveOrigFace(orFaceName.c_str());

//    //MP.SaveSharpFeatures(sharpName.c_str());

//    //tri_mesh.SaveTriMesh()
}

void DoBatchProcess ()
{
    typename MeshPrepocess<FieldTriMesh>::BatchParam BPar;
    BPar.DoRemesh=do_remesh;
    BPar.feature_erode_dilate=feature_erode_dilate;
    BPar.remesher_aspect_ratio=remesher_aspect_ratio;
    BPar.remesher_iterations=remesher_iterations;
    BPar.remesher_termination_delta=remesher_termination_delta;
    BPar.SharpFactor=SharpFactor;
    BPar.sharp_feature_thr=sharp_feature_thr;
    BPar.surf_dist_check=surf_dist_check;
    BPar.UpdateSharp=(!(has_features || has_features_fl));
    MeshPrepocess<FieldTriMesh>::BatchProcess(tri_mesh,BPar,FieldParam);
    drawfield=true;
}

void DoAutoRemesh()
{
    RemPar.iterations   = remesher_iterations;
    RemPar.targetAspect = remesher_aspect_ratio;
    RemPar.targetDeltaFN= remesher_termination_delta;
    RemPar.creaseAngle  = sharp_feature_thr;
    RemPar.erodeDilate  = feature_erode_dilate;

    //RemPar.userSelectedCreases = true;
    RemPar.surfDistCheck = surf_dist_check;

    AutoRemesher<FieldTriMesh>::RemeshAdapt(tri_mesh,RemPar);
    //AutoRemesher<FieldTriMesh>::RemeshAdapt(tri_mesh,RemPar);

    tri_mesh.InitEdgeType();
    tri_mesh.InitFeatureCoordsTable();
//    tri_mesh.InitSharpFeatures(sharp_feature_thr);
//    tri_mesh.ErodeDilate(feature_erode_dilate);

    //MeshPrepocess<FieldTriMesh>::AutoRemesh(tri_mesh,RemPar);
//    Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
}

void TW_CALL CleanMesh(void *)
 {
   MeshPrepocess<FieldTriMesh>::SolveGeometricArtifacts(tri_mesh);
 }

void TW_CALL AutoRemesh(void *)
{
    DoAutoRemesh();
}

void TW_CALL InitSharpFeatures(void *)
{
   InitSharp();
}

void TW_CALL RefineIfNeeded(void *)
{
   // tri_mesh.RefineIfNeeded();
    MeshPrepocess<FieldTriMesh>::RefineIfNeeded(tri_mesh);
    //MP.RefineIfNeeded();
    //Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
}

void TW_CALL BatchProcess(void *)
{
    DoBatchProcess();
}

void TW_CALL SmoothField(void *)
{
//    //tri_mesh.SplitFolds();
//    MeshPrepocess<FieldTriMesh>::SplitFolds(tri_mesh);
//    //tri_mesh.RemoveFolds();
//    MeshPrepocess<FieldTriMesh>::RemoveFolds(tri_mesh);
//    //tri_mesh.RemoveSmallComponents();
//    MeshPrepocess<FieldTriMesh>::RemoveSmallComponents(tri_mesh);
//    MP.SplitFolds();
//    MP.RemoveFolds();
//    MP.RemoveSmallComponents();

    //vcg::tri::io::ExporterPLY<FieldTriMesh>::Save(tri_mesh,"test0.ply");
    //tri_mesh.SmoothField(FieldParam);
    MeshFieldSmoother<FieldTriMesh>::SmoothField(tri_mesh,FieldParam);
    //vcg::tri::io::ExporterPLY<FieldTriMesh>::Save(tri_mesh,"test1.ply");
    drawfield=true;

    //Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
}



#ifdef MIQ_QUADRANGULATE
void TW_CALL MiqQuadrangulate(void *)
{

    vcg::tri::MiQParametrizer<FieldTriMesh>::MIQParametrize(tri_mesh,MiqP);//,MaxFeature);
    vcg::tri::Quadrangulator<FieldTriMesh,PMesh> Quadr;

    std::cout<<"Quadrangulating"<<std::endl;
    FieldTriMesh splitted_mesh;
    vcg::tri::Append<FieldTriMesh,FieldTriMesh>::Mesh(splitted_mesh,tri_mesh);
    std::vector< std::vector< short int> > AxisUV;

    Quadr.Quadrangulate(splitted_mesh,quad_mesh,AxisUV,false);
    quadrangulated=true;
    quad_mesh.UpdateAttributes();
    vcg::tri::Allocator<PMesh>::CompactEveryVector(quad_mesh);
    std::cout<<"Num Faces: "<<quad_mesh.face.size()<<std::endl;
    vcg::tri::io::ExporterOBJ<PMesh>::Save(quad_mesh,"quadrangulated.obj",vcg::tri::io::Mask::IOM_BITPOLYGONAL);
}
#endif



void TW_CALL SaveData(void *)
{
    SaveAllData();
}

void TW_CALL AutoSetupField(void *)
{
    MeshFieldSmoother<FieldTriMesh>::AutoSetupParam(tri_mesh,FieldParam);
}

void TW_CALL ErodeDilateFeatureStep(void *)
{
    tri_mesh.ErodeDilate(feature_erode_dilate);
    //MP.ErodeDilate(feature_erode_dilate);
}

void SetFieldBarSizePosition(QWidget *w)
{
    int params[2];
    params[0] = QTDeviceWidth(w) / 3;
    params[1] = QTDeviceHeight(w) / 1.8;
    TwSetParam(barQuad, NULL, "size", TW_PARAM_INT32, 2, params);
    params[0] = QTLogicalToDevice(w, 10);
    params[1] = 30;//QTDeviceHeight(w) - params[1] - QTLogicalToDevice(w, 10);
    TwSetParam(barQuad, NULL, "position", TW_PARAM_INT32, 2, params);
}

void InitFieldBar(QWidget *w)
{
    barQuad = TwNewBar("QuadWild");

    SetFieldBarSizePosition(w);

	TwAddVarRW(barQuad,"sharp_feature_thr",TW_TYPE_DOUBLE, &sharp_feature_thr," label='Sharp Degree'");
    TwAddVarRW(barQuad,"LimitConcave",TW_TYPE_DOUBLE, &tri_mesh.LimitConcave," label='Limit Concave'");
    //TwAddVarRW(barQuad,"LimitConcave",TW_TYPE_DOUBLE, &MP.LimitConcave," label='Limit Concave'");

    TwAddButton(barQuad,"CleanMesh",CleanMesh,0,"label='CleanMesh'");

    TwAddButton(barQuad,"SetSharp",InitSharpFeatures,0,"label='InitSharp'");

	TwAddVarRW(barQuad,"ErodeDilSteps",TW_TYPE_INT32,&feature_erode_dilate,"label='ErodeDilateSteps'");

    TwAddButton(barQuad,"Erode Dilate",ErodeDilateFeatureStep,0,"label='ErodeDilateSharp'");

    TwAddButton(barQuad,"AutoRemesh",AutoRemesh,0,"label='AutoRemesh'");

    TwAddButton(barQuad,"Refine",RefineIfNeeded,0,"label='Refine if needed'");

    TwAddVarRW(barQuad,"Alpha",TW_TYPE_DOUBLE, &FieldParam.alpha_curv," label='Alpha Curvature'");
    TwAddVarRW(barQuad,"HardCT",TW_TYPE_DOUBLE, &FieldParam.curv_thr," label='Hard Curv Thr'");
    TwAddVarRW(barQuad,"CurvRing",TW_TYPE_INT32,&FieldParam.curvRing,"label='Curvature Ring'");

    TwEnumVal smoothmodes[3] = {
        {vcg::tri::SMMiq,"MIQ"},
        {vcg::tri::SMNPoly,"NPoly"},
        {vcg::tri::SMIterative,"Ite"}
    };
    TwType smoothMode = TwDefineEnum("SmoothMode", smoothmodes, 3);
    TwAddVarRW(barQuad, "Smooth Mode", smoothMode, &FieldParam.SmoothM," label='Smooth Mode' ");


    TwAddButton(barQuad,"AutoSetup",AutoSetupField,0,"label='Auto Setup Field'");

    TwAddButton(barQuad,"ComputeField",SmoothField,0,"label='Compute Field'");

    TwAddButton(barQuad,"BatchProcess",BatchProcess,0,"label='Batch Process'");

    TwAddButton(barQuad,"SaveData",SaveData,0,"label='Save Data'");

#ifdef MIQ_QUADRANGULATE
    TwAddSeparator(barQuad,"","");
    TwAddVarRW(barQuad,"Gradient",TW_TYPE_DOUBLE, &MiqP.gradient," label='Gradient'");
    TwAddVarRW(barQuad,"Direct Round",TW_TYPE_BOOLCPP, &MiqP.directRound," label='Direct Round'");
    TwAddVarRW(barQuad,"Round Singularities",TW_TYPE_BOOLCPP, &MiqP.round_singularities," label='Round Singularities'");
    TwAddVarRW(barQuad,"IsotropyVsAlign",TW_TYPE_DOUBLE, &MiqP.miqAnisotropy," label='Isotropy Vs Align'");
    TwAddVarRW(barQuad,"Align Sharp",TW_TYPE_BOOLCPP, & MiqP.crease_as_feature," label='Align Sharp'");
    TwAddButton(barQuad,"Quadrangulate",MiqQuadrangulate,0,"label='Miq Quadrangulate'");
#endif

}

//void BatchProcess ()
//{
////    vcg::tri::Hole<FieldTriMesh>::EarCuttingFill<vcg::tri::TrivialEar<FieldTriMesh> >(tri_mesh,6);
//    bool Oriented,Orientable;
//    vcg::tri::Clean<FieldTriMesh>::OrientCoherentlyMesh(tri_mesh,Oriented,Orientable);
//    if (!Orientable)
//        std::cout<<"WARNING MESH NON ORIENTABLE"<<std::endl;

////    size_t numDup=MeshPrepocess<FieldTriMesh>::NumDuplicatedV(tri_mesh);//tri_mesh.NumDuplicatedV();
////    if (numDup>0)std::cout<<"WARNING DUPLICATED VERTICES BEFORE AUTO REMESH!"<<std::endl;

//    if (do_remesh)
//        DoAutoRemesh();

//    vcg::tri::Clean<FieldTriMesh>::OrientCoherentlyMesh(tri_mesh,Oriented,Orientable);
//    if (!Orientable)
//        std::cout<<"WARNING MESH NON ORIENTABLE"<<std::endl;

//    for (size_t i=0;i<tri_mesh.face.size();i++)
//        tri_mesh.face[i].IndexOriginal=i;

//    //tri_mesh.SolveGeometricIssues();
//    MeshPrepocess<FieldTriMesh>::SolveGeometricIssues(tri_mesh);
//    //MP.SolveGeometricIssues();


//	std::string projM=pathM;
//    size_t indexExt=projM.find_last_of(".");
//    projM=projM.substr(0,indexExt);
//    std::string meshName=projM+std::string("_remeshed.obj");
//    std::cout<<"Saving remeshed Mesh TO:"<<meshName.c_str()<<std::endl;
//    tri_mesh.SaveTriMesh(meshName.c_str());

//    vcg::tri::Allocator<FieldTriMesh>::CompactEveryVector(tri_mesh);

//    if (!(has_features || has_features_fl))
//    {
//        InitSharp();
//        tri_mesh.ErodeDilate(feature_erode_dilate);
//        //MP.ErodeDilate(feature_erode_dilate);
//    }

//    //tri_mesh.RefineIfNeeded();
//    MeshPrepocess<FieldTriMesh>::RefineIfNeeded(tri_mesh);
//    //tri_mesh.SolveGeometricIssues();
//    MeshPrepocess<FieldTriMesh>::SolveGeometricIssues(tri_mesh);
//    //MP.RefineIfNeeded();
//    //MP.SolveGeometricIssues();

//    //FIELD SMOOTH
//    bool UseNPoly=tri_mesh.SufficientFeatures(SharpFactor);
//    ScalarType SharpL=tri_mesh.SharpLenght();
//    std::cout<<"Sharp Lenght 0: "<<SharpL<<std::endl;
//    if (UseNPoly)
//    {
//        std::cout<<"Using NPoly"<<std::endl;
//        FieldParam.SmoothM=SMNPoly;
//    }
//    else
//    {
//        std::cout<<"Using Comiso"<<std::endl;
//        FieldParam.SmoothM=SMMiq;
//        FieldParam.alpha_curv=0.3;
//    }

//    std::string projM=pathM;
//    size_t indexExt=projM.find_last_of(".");
//    projM=projM.substr(0,indexExt);
//    std::string meshName=projM+std::string("_rem.obj");
//    std::string sharpName=projM+std::string("_rem.sharp");
//    std::cout<<"Saving Mesh TO:"<<meshName.c_str()<<std::endl;
//    std::cout<<"Saving Sharp TO:"<<sharpName.c_str()<<std::endl;

//    std::cout << "[fieldComputation] Smooth Field Computation..." << std::endl;
//    //tri_mesh.SplitFolds();
//    MeshPrepocess<FieldTriMesh>::SplitFolds(tri_mesh);

//    SharpL=tri_mesh.SharpLenght();
//    //SharpL=MP.SharpLenght();
//    std::cout<<"Sharp Lenght 0.0: "<<SharpL<<std::endl;
//    //tri_mesh.RemoveFolds();
//    MeshPrepocess<FieldTriMesh>::RemoveFolds(tri_mesh);

//    SharpL=tri_mesh.SharpLenght();
//    //SharpL=MP.SharpLenght();
//    std::cout<<"Sharp Lenght 0.1: "<<SharpL<<std::endl;
//    //tri_mesh.SolveGeometricIssues();
//    MeshPrepocess<FieldTriMesh>::SolveGeometricIssues(tri_mesh);
//    //MP.SolveGeometricIssues();

//    SharpL=tri_mesh.SharpLenght();
//    //SharpL=MP.SharpLenght();
//    std::cout<<"Sharp Lenght 0.2: "<<SharpL<<std::endl;
//    //tri_mesh.RemoveSmallComponents();
//    MeshPrepocess<FieldTriMesh>::RemoveSmallComponents(tri_mesh);
//    //MP.RemoveSmallComponents();
//    SharpL=tri_mesh.SharpLenght();
//    //SharpL=MP.SharpLenght();
//    std::cout<<"Sharp Lenght 0.3: "<<SharpL<<std::endl;

//    SharpL=tri_mesh.SharpLenght();
//    //SharpL=MP.SharpLenght();
//    std::cout<<"Sharp Lenght 1: "<<SharpL<<std::endl;
//    //tri_mesh.SmoothField(FieldParam);
//    MeshFieldSmoother<FieldTriMesh>::SmoothField(tri_mesh,FieldParam);

//    SharpL=tri_mesh.SharpLenght();
//    //SharpL=MP.SharpLenght();
//    std::cout<<"Sharp Lenght 2: "<<SharpL<<std::endl;

//    std::string fieldName=projM+std::string("_rem.rosy");
//    std::cout<<"Saving Field TO:"<<fieldName.c_str()<<std::endl;

//    tri_mesh.SaveTriMesh(meshName.c_str());
//    tri_mesh.SaveSharpFeatures(sharpName.c_str());
//    //MP.SaveSharpFeatures(sharpName.c_str());
//    tri_mesh.SaveField(fieldName.c_str());

//    std::string orFaceName=projM+std::string("_rem_origf.txt");
//    std::cout<<"Saving Field TO:"<<orFaceName.c_str()<<std::endl;
//    tri_mesh.SaveOrigFace(orFaceName.c_str());
//    drawfield=true;
////    //SAVE
//    SaveAllData();
//}




GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
    hasToPick=false;
    bool AllQuad=false;
    bool Loaded=tri_mesh.LoadTriMesh(pathM,AllQuad);
    FieldParam.alpha_curv=alpha;
    FieldParam.curv_thr=0.8;

    if (!Loaded)
    {
        std::cout<<"Error Loading Mesh"<<std::endl;
        exit(0);
    }
    if (AllQuad)
    {
        drawfield=true;
    }
    std::cout<<"Loaded "<<tri_mesh.face.size()<<" faces "<<std::endl;
    std::cout<<"Loaded "<<tri_mesh.vert.size()<<" vertices "<<std::endl;

    glWrap.m=&tri_mesh;

    tri_mesh.UpdateDataStructures();

    tri_mesh.LimitConcave=0;
    //MP.LimitConcave=0;
    if (has_features)
    {
        std::cout<<"*** Loading SHARP FEATURES ***"<<std::endl;
        bool HasRead=tri_mesh.LoadSharpFeatures(pathS);
        //bool HasRead=MP.LoadSharpFeatures(pathS);
        assert(HasRead);
    }

    if (has_features_fl)
    {
        std::cout<<"*** Loading SHARP FEATURES FL***"<<std::endl;
        bool HasRead=tri_mesh.LoadSharpFeaturesFL(pathS);
        //bool HasRead=MP.LoadSharpFeaturesFL(pathS);
        assert(HasRead);
    }

//    if ((do_batch)&&(!has_features))
//    {
////        bool SufficientFeatures=mesh.SufficientFeatures(SharpFactor);
////        if (SufficientFeatures)

//        DoBatchProcess();
//        SaveAllData();
//        exit(0);
//    }
//    if ((do_batch)&&(has_features))
//    {
//        //tri_mesh.PrintSharpInfo();
//        DoBatchProcess();
//        SaveAllData();
//        exit(0);
//    }
    if (do_batch)
    {

        DoBatchProcess();
        SaveAllData();
        exit(0);
    }

    //remeshed_mesh.UpdateDataStructures();
    //Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
}



void GLWidget::initializeGL ()
{
    //initialize Glew
    glewInit();
    //CaptInt.GLInit( GLWidget::width(),GLWidget::height());
    glClearColor(0, 0, 0, 0);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
}


void GLWidget::resizeGL (int w, int h)
{
    glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    TwWindowSize(w, h);
    InitFieldBar(this);
    initializeGL();
}

void GLWidget::paintGL ()
{

    //    if (RType!=OldRType)
    //    {
    //        PatchDeco.ColorPatches(RType);
    //        OldRType=RType;
    //    }

    glClearColor(255,255,255,255);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, GLWidget::width()/(float)GLWidget::height(), 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,3.5f,   0,0,0,   0,1,0);
    track.center=vcg::Point3f(0, 0, 0);
    track.radius= 1;
    track.GetView();

    glPushMatrix();
    track.Apply();
    glPushMatrix();

    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    if(tri_mesh.vert.size()>0)
    {
        vcg::glScale(2.0f/tri_mesh.bbox.Diag());
        glTranslate(-tri_mesh.bbox.Center());
        //tri_mesh.GLDrawSharpEdges();
        GLDrawSharpEdges(tri_mesh);
        //MP.GLDrawSharpEdges();
        glWrap.Draw(vcg::GLW::DMFlatWire,vcg::GLW::CMPerFace,vcg::GLW::TMNone);
        //glWrap.Draw(vcg::GLW::DMSmooth,vcg::GLW::CMNone,vcg::GLW::TMNone);
    }

#ifdef MIQ_QUADRANGULATE
    if (quadrangulated)
    {
        quad_mesh.GLDraw();
    }
#endif

    if (drawfield)
    {
        vcg::GLField<FieldTriMesh>::GLDrawFaceField(tri_mesh,false,false,0.007);
        vcg::GLField<FieldTriMesh>::GLDrawSingularity(tri_mesh);
    }

//    if(hasToPick)
//    {
//        hasToPick=false;
//        typename FieldTriMesh::CoordType pp;
//        if(vcg::Pick<typename FieldTriMesh::CoordType>(xMouse,yMouse,pp))
//        {
//            typename FieldTriMesh::CoordType closPt,bary;
//            typename FieldTriMesh::ScalarType minD;
//            typename FieldTriMesh::FaceType *f=vcg::tri::GetClosestFaceBase(tri_mesh,Gr,pp,tri_mesh.bbox.Diag(),minD,closPt);
//            vcg::InterpolationParameters(*f,closPt,bary);
//            size_t EdgeI=1;
//            if ((bary.Y()<bary.X())&&(bary.Y()<bary.Z()))EdgeI=2;
//            if ((bary.Z()<bary.X())&&(bary.Z()<bary.Y()))EdgeI=0;

////            FieldTriMesh::FaceType *fOpp=f->FFp(EdgeI);
////            int eOpp=f->FFi(EdgeI);

//            if (f->IsFaceEdgeS(EdgeI))
//            {
//                tri_mesh.ClearSharp((*f),EdgeI);
//                //MP.ClearSharp((*f),EdgeI);
////                {
////                f->ClearFaceEdgeS(EdgeI);
////                if (fOpp!=f)
////                    fOpp->ClearFaceEdgeS(eOpp);
//            }else
//            {

//                tri_mesh.SetSharp((*f),EdgeI);
//                //MP.SetSharp((*f),EdgeI);
////                f->SetFaceEdgeS(EdgeI);
////                if (fOpp!=f)
////                    fOpp->SetFaceEdgeS(eOpp);
//            }
//        }
//    }

    glPopMatrix();
    glPopMatrix();


    TwDraw();

}

void GLWidget::keyReleaseEvent (QKeyEvent * e)
{
    e->ignore ();
    if (e->key () == Qt::Key_Control)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ControlModifier));
    if (e->key () == Qt::Key_Shift)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
    if (e->key () == Qt::Key_Alt) track.ButtonUp (QT2VCG (Qt::NoButton, Qt::AltModifier));
    updateGL ();
}


void GLWidget::keyPressEvent (QKeyEvent * e)
{
    e->ignore ();
    if (e->key () == Qt::Key_Control) track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ControlModifier));
    if (e->key () == Qt::Key_Shift)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
    if (e->key () == Qt::Key_Alt)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::AltModifier));

    TwKeyPressQt(e);
    updateGL ();
}

void GLWidget::mousePressEvent (QMouseEvent * e)
{
    if(!TwMousePressQt(this,e))
    {
        e->accept ();
        setFocus ();
        track.MouseDown(QT2VCG_X(this, e), QT2VCG_Y(this, e), QT2VCG (e->button (), e->modifiers ()));
    }
    updateGL ();
}

void GLWidget::mouseMoveEvent (QMouseEvent * e)
{
    if (e->buttons ()) {
        track.MouseMove(QT2VCG_X(this, e), QT2VCG_Y(this, e));
        updateGL ();
    }
    TwMouseMotion(QTLogicalToDevice(this, e->x()), QTLogicalToDevice(this, e->y()));
}

void GLWidget::mouseDoubleClickEvent (QMouseEvent * e)
{
    if (e->buttons ())
    {
        xMouse=QT2VCG_X(this, e);
        yMouse=QT2VCG_Y(this, e);
        //pointToPick=Point2i(e->x(),height()-e->y());
        //pointToPick=Point2i(xMouse,yMouse);
        hasToPick=true;
        updateGL ();
    }
    updateGL();
}

void GLWidget::mouseReleaseEvent (QMouseEvent * e)
{
    track.MouseUp(QT2VCG_X(this, e), QT2VCG_Y(this, e), QT2VCG(e->button (), e->modifiers ()));
    TwMouseReleaseQt(this,e);
    updateGL ();
}

void GLWidget::wheelEvent (QWheelEvent * e)
{
    const int WHEEL_STEP = 120;
    track.MouseWheel (e->delta () / float (WHEEL_STEP), QTWheel2VCG (e->modifiers ()));
    updateGL ();
}
