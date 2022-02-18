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

#include "glwidget.h"
#include <wrap/qt/trackball.h>
#include <wrap/gl/picking.h>
#include <wrap/gl/trimesh.h>
#include <wrap/qt/anttweakbarMapper.h>
#include <wrap/gl/gl_field.h>
#include <QDir>
#include "tracing/GL_vert_field_graph.h"
#include "tracing/patch_tracer.h"
#include "tracing/tracer_interface.h"
#include "tracing/GL_mesh_drawing.h"
#include "tracing/GL_metamesh.h"
//#include <vcg/complex/algorithms/parametrization/local_para_smooth.h>

TwBar *bar;
char * filename;/// filename of the mesh to load
TraceMesh mesh,problematic_mesh;     /// the active mesh instance
vcg::GlTrimesh<TraceMesh> glWrap;    /// the active mesh opengl wrapper
vcg::GlTrimesh<TraceMesh> glWrapProblem;
vcg::Trackball track;     /// the active manipulator
GLW::DrawMode drawmode=GLW::DMFlatWire;     /// the current drawmode

std::string pathM,pathF,pathS,pathOF,pathProject;

//bool to_pick=false;
//int xMouse,yMouse;

bool drawSharpF=true;
bool drawSing=true;
bool drawField=true;
bool drawPaths=true;
bool drawPathNodes=false;
bool drawTwins=false;

bool drawConcaveEmitters=false;
bool drawConcaveReceivers=false;
bool drawFlatEmitters=false;
bool drawFlatReceivers=false;
bool drawNarrowEmitters=false;
bool drawNarrowReceivers=false;
bool drawInvalidated=false;
bool drawUnsatisfied=false;
bool drawCorners=false;
bool drawChoosenEmitters=false;
bool drawChoosenReceivers=false;
bool drawMetaMesh=false;
bool drawNarrowCandidates=false;

//bool splitted=false;
bool save_setup=false;
bool save_csv=false;

bool batch_process=false;
bool has_features=false;
bool has_original_faces=false;
bool add_only_needed=false;
bool meta_mesh_collapse=true;
//bool interleave_removal=true;
//bool interleave_smooth=false;
bool final_removal=true;
bool force_split=false;
bool drawProblematicsOnly=false;
int CurrNum=0;
int VisTraces=-1;

bool subdivide_when_save=false;

std::vector<size_t> ConcaveEmittersNode,ConcaveReceiversNode,
FlatEmittersNode,FlatReceiversNode,
NarrowEmittersNode,NarrowReceiversNode,
ProblematicNodes,UnsatisfiedNodes,
ChoosenEmittersNode,ChoosenReceiversNode,
TraceableFlatNode;

VertexFieldGraph<TraceMesh> VGraph(mesh);
GLVertGraph<TraceMesh> GLGraph(VGraph);

typedef PatchTracer<TraceMesh> TracerType;
TracerType PTr(VGraph);

std::vector<std::vector<size_t> > CurrCandidates;
std::vector<bool> ChosenIsLoop;
std::vector<std::vector<size_t> > ChosenCandidates;
std::vector<bool> DiscardedIsLoop;
std::vector<std::vector<size_t> > DiscardedCandidates;
std::vector<typename TraceMesh::CoordType> PatchCornerPos;
TraceMesh::ScalarType Drift=100;

enum PatchColorMode{CMPatchNone, CMPatchCol, CMPatchValence,
                    CMPatchTopo,CMPatchLengthDist,CMPatchLenghtVar,
                    CMArapDistortion,CMCClarkability,CMExpValence};

PatchColorMode CurrPatchMode=CMPatchNone;
PatchColorMode OldPatchMode=CMPatchNone;

//float BatchSample=-1;
//float BatchDrift=-1;
//int BatchSplit=-1;
//int BatchIncreaseValRem=-1;
//int BatchIrregularRem=-1;
//float BatchDistortionL=-1;


void SaveSetupFile(const std::string pathProject,
                   const size_t CurrNum)
{
    std::string pathSetupFinal=pathProject;
    pathSetupFinal=pathSetupFinal+"_p"+std::to_string(CurrNum)+".setup";

    FILE *f=fopen(pathSetupFinal.c_str(),"wt");
    assert(f!=NULL);

    //    fprintf(f,"Srate %f\n",PTr.sample_ratio);
    //    fprintf(f,"Drift %f\n",Drift);

    //    if (PTr.split_on_removal)
    //        fprintf(f,"Split 1\n");
    //    else
    //        fprintf(f,"Split 0\n");


    //    FILE *f=fopen(path.c_str(),"rt");
    //    assert(f!=NULL);

    //    float Driftf;
    //    fscanf(f,"Drift %f\n",&Driftf);
    //    Drift=(TraceMesh::ScalarType)Driftf;
    //    std::cout<<"DRIFT "<<Drift<<std::endl;
    fprintf(f,"Drift %f\n",Drift);

    //    float SRatef;
    //    fscanf(f,"Srate %f\n",&SRatef);
    //    PTr.sample_ratio=(TraceMesh::ScalarType)SRatef;
    //    std::cout<<"SAMPLE RATE "<<PTr.sample_ratio<<std::endl;
    fprintf(f,"Srate %f\n",PTr.sample_ratio);

    //    int IntVar=0;
    //    fscanf(f,"SplitOnRem %d\n",&IntVar);
    //    std::cout<<"SPLIT "<<IntVar<<std::endl;
    //    if (IntVar==0)
    //        PTr.split_on_removal=false;
    //    else
    //        PTr.split_on_removal=true;

    if (PTr.split_on_removal)
        fprintf(f,"SplitOnRem 1\n");
    else
        fprintf(f,"SplitOnRem 0\n");

    //fscanf(f,"MaxVal %d\n",&PTr.MaxVal);
    fprintf(f,"MaxVal %d\n",PTr.MaxVal);
    //std::cout<<"INCREASE VAL "<<PTr.MaxVal<<std::endl;

    //    float CCbility;
    //    fscanf(f,"CCability %f\n",&CCbility);
    //    PTr.CClarkability=(TraceMesh::ScalarType)CCbility;
    //    std::cout<<"CCABILITY "<<PTr.CClarkability<<std::endl;
    fprintf(f,"CCability %f\n",PTr.CClarkability);

    //    fscanf(f,"MatchVal %d\n",&IntVar);
    //    std::cout<<"MATCH VAL "<<IntVar<<std::endl;
    //    if (IntVar==0)
    //        PTr.match_valence=false;
    //    else
    //        PTr.match_valence=true;
    if (PTr.match_valence)
        fprintf(f,"MatchVal 1\n");
    else
        fprintf(f,"MatchVal 0\n");

    //    fscanf(f,"AddNeed %d\n",&IntVar);
    //    std::cout<<"ADD NEED "<<IntVar<<std::endl;
    //    if (IntVar==0)
    //        add_only_needed=false;
    //    else
    //        add_only_needed=true;
    if (add_only_needed)
        fprintf(f,"AddNeed 1\n");
    else
        fprintf(f,"AddNeed 0\n");

    //    fscanf(f,"InterleaveRem %d\n",&IntVar);
    //    std::cout<<"INTERLEAVE REMOVE "<<IntVar<<std::endl;
    //    if (IntVar==0)
    //        interleave_removal=false;
    //    else
    //        interleave_removal=true;
    //    if (interleave_removal)
    //        fprintf(f,"InterleaveRem 1\n");
    //    else
    //        fprintf(f,"InterleaveRem 0\n");

    //    fscanf(f,"InterleaveSmooth %d\n",&IntVar);
    //    std::cout<<"INTERLEAVE SMOOTH "<<IntVar<<std::endl;
    //    if (IntVar==0)
    //        interleave_smooth=false;
    //    else
    //        interleave_smooth=true;
    //    if (interleave_smooth)
    //        fprintf(f,"InterleaveSmooth 1\n");
    //    else
    //        fprintf(f,"InterleaveSmooth 0\n");

    //    fscanf(f,"FinalRem %d\n",&IntVar);
    //    std::cout<<"FINAL REMOVE "<<IntVar<<std::endl;
    //    if (IntVar==0)
    //        final_removal=false;
    //    else
    //        final_removal=true;
    if (final_removal)
        fprintf(f,"FinalRem 1\n");
    else
        fprintf(f,"FinalRem 0\n");


    if (force_split)
        fprintf(f,"Force Split 1\n");
    else
        fprintf(f,"Force Split 0\n");

    if (subdivide_when_save)
        fprintf(f,"Subd 1\n");
    else
        fprintf(f,"Subd 0\n");

    if (meta_mesh_collapse)
        fprintf(f,"MetaCollapse 1\n");
    else
        fprintf(f,"MetaCollapse 0\n");

    //    if (PTr.avoid_increase_valence)
    //        fprintf(f,"IncreaseValRem 1\n");
    //    else
    //        fprintf(f,"IncreaseValRem 0\n");

    //    if (PTr.avoid_collapse_irregular)
    //        fprintf(f,"IrregularRem 1\n");
    //    else
    //        fprintf(f,"IrregularRem 0\n");

    //fprintf(f,"DistortionL %f\n",PTr.max_lenght_distortion);

    fclose(f);
}





void LoadSetupFile(std::string path)
{
    FILE *f=fopen(path.c_str(),"rt");
    assert(f!=NULL);

    float Driftf;
    fscanf(f,"Drift %f\n",&Driftf);
    Drift=(TraceMesh::ScalarType)Driftf;
    std::cout<<"DRIFT "<<Drift<<std::endl;

    float SRatef;
    fscanf(f,"Srate %f\n",&SRatef);
    PTr.sample_ratio=(TraceMesh::ScalarType)SRatef;
    std::cout<<"SAMPLE RATE "<<PTr.sample_ratio<<std::endl;

    int IntVar=0;
    fscanf(f,"SplitOnRem %d\n",&IntVar);
    std::cout<<"SPLIT "<<IntVar<<std::endl;
    if (IntVar==0)
        PTr.split_on_removal=false;
    else
        PTr.split_on_removal=true;

    fscanf(f,"MaxVal %d\n",&PTr.MaxVal);
    std::cout<<"INCREASE VAL "<<PTr.MaxVal<<std::endl;

    float CCbility;
    fscanf(f,"CCability %f\n",&CCbility);
    PTr.CClarkability=(TraceMesh::ScalarType)CCbility;
    std::cout<<"CCABILITY "<<PTr.CClarkability<<std::endl;

    fscanf(f,"MatchVal %d\n",&IntVar);
    std::cout<<"MATCH VAL "<<IntVar<<std::endl;
    if (IntVar==0)
        PTr.match_valence=false;
    else
        PTr.match_valence=true;

    fscanf(f,"AddNeed %d\n",&IntVar);
    std::cout<<"ADD NEED "<<IntVar<<std::endl;
    if (IntVar==0)
        add_only_needed=false;
    else
        add_only_needed=true;

    //    fscanf(f,"InterleaveRem %d\n",&IntVar);
    //    std::cout<<"INTERLEAVE REMOVE "<<IntVar<<std::endl;
    //    if (IntVar==0)
    //        interleave_removal=false;
    //    else
    //        interleave_removal=true;

    //    fscanf(f,"InterleaveSmooth %d\n",&IntVar);
    //    std::cout<<"INTERLEAVE SMOOTH "<<IntVar<<std::endl;
    //    if (IntVar==0)
    //        interleave_smooth=false;
    //    else
    //        interleave_smooth=true;

    fscanf(f,"FinalRem %d\n",&IntVar);
    std::cout<<"FINAL REMOVE "<<IntVar<<std::endl;
    if (IntVar==0)
        final_removal=false;
    else
        final_removal=true;


    fscanf(f,"ForceSplit %d\n",&IntVar);
    std::cout<<"FORCE SPLIT "<<IntVar<<std::endl;
    if (IntVar==0)
        force_split=false;
    else
        force_split=true;

    fscanf(f,"Subd %d\n",&IntVar);
    std::cout<<"SUBDIVIDE WHEN SAVE "<<IntVar<<std::endl;
    if (IntVar==0)
        subdivide_when_save=false;
    else
        subdivide_when_save=true;

    fscanf(f,"MetaCollapse %d\n",&IntVar);
    std::cout<<"META COLLAPSE "<<IntVar<<std::endl;
    if (IntVar==0)
        meta_mesh_collapse=false;
    else
        meta_mesh_collapse=true;

    //    if ((batch_process)&&(BatchSample>0))
    //        PTr.sample_ratio=BatchSample;

    //    if ((batch_process)&&(BatchDrift>0))
    //        Drift=BatchDrift;

    //    if ((batch_process)&&(BatchSplit==0))
    //        PTr.split_on_removal=false;

    //    if ((batch_process)&&(BatchSplit==1))
    //        PTr.split_on_removal=true;
}

void FindCurrentNum()
{
    std::string BasePath=pathProject;
    BasePath.resize(BasePath.find_last_of("/")+1);
    BasePath="./"+BasePath;
    //std::cout<<"basepath "<<BasePath.c_str()<<std::endl;
    QDir directory(BasePath.c_str());

    QFile f(pathM.c_str());
    QFileInfo fileInfo(f.fileName());
    QString filename(fileInfo.fileName());
    std::string NameFile=filename.toStdString();
    //std::cout<<"namefile "<<NameFile.c_str()<<std::endl;
    NameFile.resize(NameFile.find_last_of("."));
    std::string Filter=NameFile+"_p*.obj";
    //std::cout<<"filter "<<Filter.c_str()<<std::endl;

    //    TestPath+="*.obj";
    QStringList projectFiles = directory.entryList(QStringList() <<Filter.c_str(),QDir::Files);
    CurrNum=0;

    foreach(QString filename, projectFiles) {
        int Num;
        std::string TestScan=NameFile+"_p%d.obj";
        sscanf (filename.toStdString().c_str(),TestScan.c_str(),&Num);
        CurrNum=std::max(CurrNum,(Num+1));
        //        std::cout<<"test "<<Num<<std::endl;
        //        std::cout<<"test "<<filename.toStdString().c_str()<<std::endl;
        //do whatever you need to do
    }
}


void UpdatePatchColor()
{
    if (CurrPatchMode==CMPatchNone)
        vcg::tri::UpdateColor<TraceMesh>::PerFaceConstant(mesh,vcg::Color4b(192,192,192,255));
    if (CurrPatchMode==CMPatchCol)
        PTr.ColorByPartitions();
    //PTr.ColorByPatchQuality();
    if (CurrPatchMode==CMPatchValence)
        PTr.ColorByValence();
    if (CurrPatchMode==CMPatchTopo)
        PTr.ColorByTopology();
    if (CurrPatchMode==CMPatchLenghtVar)
        PTr.ColorByLenghtVariance();
    if (CurrPatchMode==CMPatchLengthDist)
        PTr.ColorByLenghtDistortion();
    if (CurrPatchMode==CMArapDistortion)
        PTr.ColorByArapDistortion();
    if (CurrPatchMode==CMCClarkability)
        PTr.ColorByCClarkability();
    if (CurrPatchMode==CMExpValence)
        PTr.ColorByExpValence();
}

void UpdateVisualNodes()
{

    CurrCandidates.clear();
    PTr.GetCurrCandidates(CurrCandidates);

    ChosenCandidates.clear();
    PTr.GetCurrChosen(ChosenCandidates);

    DiscardedCandidates.clear();
    PTr.GetCurrDiscarded(DiscardedCandidates);

    ChosenIsLoop.clear();
    PTr.GetCurrChosenIsLoop(ChosenIsLoop);

    DiscardedIsLoop.clear();
    PTr.GetCurrDiscardedIsLoop(DiscardedIsLoop);
    //return;

    //PTr.GetConcaveEmitters(ConcaveEmittersNode);
    PTr.GetActiveEmittersType(TVConcave,ConcaveEmittersNode);

    //PTr.GetConcaveReceivers(ConcaveReceiversNode);
    PTr.GetActiveReceiversType(TVConcave,ConcaveReceiversNode);

    //PTr.GetFlatEmitters(FlatEmittersNode);
    PTr.GetActiveEmittersType(TVFlat,FlatEmittersNode);

    //PTr.GetFlatReceivers(FlatReceiversNode);
    PTr.GetActiveReceiversType(TVFlat,FlatReceiversNode);

    //PTr.GetNarrowActiveEmitters(NarrowEmittersNode);
    PTr.GetActiveEmittersType(TVNarrow,NarrowEmittersNode);

    //PTr.GetNarrowActiveReceivers(NarrowReceiversNode);
    PTr.GetActiveReceiversType(TVNarrow,NarrowReceiversNode);

    //PTr.GetChoosenEmitters(ChoosenEmittersNode);
    PTr.GetActiveEmittersType(TVChoosen,ChoosenEmittersNode);

    //PTr.GetChoosenReceivers(ChoosenReceiversNode);
    PTr.GetActiveReceiversType(TVChoosen,ChoosenReceiversNode);

    PTr.GetUnsatisfied(UnsatisfiedNodes);

    PTr.GetVisualCornersPos(PatchCornerPos);

    PTr.GetTraceableFlatNodes(TraceableFlatNode);

    PTr.GetUnsolvedMesh(problematic_mesh);
}

void InitStructures()
{
    ConcaveEmittersNode.clear();
    ConcaveReceiversNode.clear();
    FlatEmittersNode.clear();
    FlatReceiversNode.clear();
    NarrowEmittersNode.clear();
    NarrowReceiversNode.clear();
    UnsatisfiedNodes.clear();

    //    std::string origfaceP=pathProject;
    //    origfaceP=featureCorners+"_p"+std::to_string(CurrNum)+"_origf.txt";
    //    mesh.SaveOrigFace("test0face.txt");
    PreProcessMesh(mesh);
    //    mesh.SaveOrigFace("test1face.txt");
    //    exit(0);
    VGraph.InitGraph(false);//SingPos);

    GLGraph.InitDisplacedPos();

    PTr.InitTracer(Drift,false);

    //std::cout<<"Here"<<std::endl;
    UpdateVisualNodes();
}

//void TW_CALL SmoothParam(void *)
//{
//    vcg::tri::Local_Param_Smooth<TraceMesh>::Smooth(mesh);
//}

void TW_CALL InitGraph(void *)
{
    //    vcg::tri::Local_Param_Smooth<TraceMesh>::SmoothStep(mesh,0.5);
    InitStructures();

    drawField=false;
    drawSharpF=false;
    drawSing=false;
}

//void TW_CALL TestInit(void *)
//{
//    PTr.sample_ratio=0.001;
//    PTr.CClarkability=-1;
//    //PTr.
//    PatchGeneralParameters::MinSamples()=1;
//    InitStructures();
//    PTr.DebugMsg=true;
//    PTr.TraceLoops(false);
//    PTr.UpdatePartitionsFromChoosen();
//    PTr.ColorByPartitions();
//    PTr.SmoothPatches();
//    UpdateVisualNodes();
//    drawField=false;
//    drawSharpF=true;
//    drawSing=true;
//}

void TW_CALL JoinNarrow(void *)
{
    PTr.JoinNarrowStep();
    PTr.UpdatePartitionsFromChoosen();
    PTr.ColorByPartitions();
    UpdateVisualNodes();
}

void TW_CALL JoinConcave(void *)
{
    PTr.JoinConcaveStep();
    PTr.UpdatePartitionsFromChoosen();
    PTr.ColorByPartitions();
    UpdateVisualNodes();
}

void TW_CALL AddLoops(void *)
{
    PTr.TraceLoops(false);
    PTr.UpdatePartitionsFromChoosen();
    PTr.ColorByPartitions();
    UpdateVisualNodes();
}

void TW_CALL TraceBorder(void *)
{
    PTr.JoinBoundaries(false);
    PTr.UpdatePartitionsFromChoosen();
    PTr.ColorByPartitions();
    UpdateVisualNodes();
}

void TW_CALL SmoothPathes(void *)
{
    PTr.SmoothPatches();
    PTr.GetVisualCornersPos(PatchCornerPos);
}

//void TW_CALL RemoveIteration(void *)
//{


////    PTr.SetAllRemovable();
////    PTr.match_valence=false;
////    PTr.AllowDarts=false;
////    PTr.DebugMsg=false;
////    PTr.AllowSelfGluedPatch=true;
////    PTr.MinVal=0;
////    PTr.CheckQuadrangulationLimits=false;
////    //PTr.SingleRemoveStep(true);
////    while(PTr.SingleRemoveStep(false)){}

//    //PTr.MergeContiguousPaths();
//    //PTr.ColorByPartitions();

//    InitStructures();

//    //RecursiveProcessForTexturing(PTr,Drift, add_only_needed,final_removal,true,meta_mesh_collapse,force_split);

//    RecursiveProcessForTexturingWithDarts(PTr,Drift, add_only_needed,final_removal,true,meta_mesh_collapse,force_split);


////    PTr.SetAllRemovable();
////    PTr.AllowDarts=true;
////    PTr.AllowSelfGluedPatch=true;
////    PTr.MinVal=0;
////    PTr.split_on_removal=true;
////    PTr.match_valence=false;
////    PTr.CheckQuadrangulationLimits=false;

////    //PTr.SplitIntoIntervals();

////    PTr.BatchRemovalOnMesh(true);

//    PTr.MergeContiguousPaths();

//    PTr.ColorByPartitions();

//    CurrPatchMode=CMPatchCol;
//    drawField=false;
//    drawSharpF=false;
//    drawSing=false;
//    UpdateVisualNodes();
//}

void TW_CALL BatchProcess(void *)
{
    InitStructures();

    //PTr.BatchProcess();
    PTr.BatchAddLoops(false,false,false,false);//,false,true);//,false);
    PTr.UpdatePartitionsFromChoosen();
    PTr.ColorByPartitions();
    CurrPatchMode=CMPatchCol;

    drawField=false;
    drawSharpF=false;
    drawSing=false;

    UpdateVisualNodes();
}

void TW_CALL RecursiveProcess(void *)
{
    InitStructures();
    //RecursiveProcess<TraceMesh>(PTr,Drift);
    RecursiveProcess<TracerType>(PTr,Drift, add_only_needed,final_removal,true,meta_mesh_collapse,force_split);//,interleave_smooth);
    CurrPatchMode=CMPatchCol;
    drawField=false;
    drawSharpF=false;
    drawSing=false;
    UpdateVisualNodes();
}

void TW_CALL ParametrizePatches(void *)
{
    PTr.ComputePatchesUV();
}

void TW_CALL SubdividePatches(void *)
{
    PTr.SubdivideIrrPatches();
    PTr.GetVisualCornersPos(PatchCornerPos);
    CurrPatchMode=CMPatchCol;

    CurrCandidates.clear();
    ChosenCandidates.clear();
    ChosenIsLoop.clear();
}

void TW_CALL BatchRemoval(void *)
{
    PTr.SetAllRemovable();

    if (meta_mesh_collapse)
        PTr.BatchRemovalMetaMesh();
    else
        PTr.BatchRemovalOnMesh();

    PTr.UpdatePartitionsFromChoosen(true);
    PTr.ColorByPartitions();

    CurrPatchMode=CMPatchCol;
    drawField=false;
    drawSharpF=false;
    drawSing=false;
    UpdateVisualNodes();
    PTr.GetVisualCornersPos(PatchCornerPos);
}

void LoadAll()
{
    printf("Loading the mesh \n");
    bool loadedMesh=mesh.LoadMesh(pathM);
    mesh.UpdateAttributes();
    if (!loadedMesh)
    {
        std::cout<<"*** ERROR LOADING MESH ***"<<std::endl;
        exit(0);
    }
    std::cout<<"Loaded "<<mesh.fn<<" faces and "<<mesh.vn<<" edges"<<std::endl;

    //FIELD LOAD
    bool loadedField=mesh.LoadField(pathF);
    if (!loadedField){
        std::cout<<"*** ERROR LOADING FIELD ***"<<std::endl;
        exit(0);
    }

    if (has_original_faces)
    {
        bool loadedOrigF=mesh.LoadOrigFaces(pathOF);

        if (!loadedOrigF){
            std::cout<<"*** ERROR LOADING ORIGINAL FACES ***"<<std::endl;
            exit(0);
        }
    }

    if (has_features){
        bool loadedFeatures=mesh.LoadSharpFeatures(pathS);
        if (!loadedFeatures){
            std::cout<<"*** ERROR LOADING FEATURES ***"<<std::endl;
            exit(0);
        }
    }

    mesh.SolveGeometricIssues();
    mesh.UpdateSharpFeaturesFromSelection();
}

void Reload()
{
    drawSing=true;
    drawField=true;
    drawPaths=false;
    LoadAll();
    CurrPatchMode=CMPatchNone;
}

void TW_CALL ReloadAll(void *)
{
    Reload();
}


std::vector<size_t> Emitter,Receiver,Disabled;

bool testDrawEmitter=false;
bool testdrawReceiver=false;
bool testdrawDisables=false;


void TW_CALL TestNarrowNarrow(void *)
{
    PTr.TestGetNodes(TVNarrow,
                     TVNarrow,
                     TraceDirect,
                     Emitter,Receiver,Disabled);
}

void TW_CALL TestNarrowConcave(void *)
{
    PTr.TestGetNodes(TVNarrow,
                     TVConcave,
                     TraceDirect,
                     Emitter,Receiver,Disabled);
}

void TW_CALL TestNarrowFlat(void *)
{
    PTr.TestGetNodes(TVNarrow,
                     TVFlat,
                     TraceDirect,
                     Emitter,Receiver,Disabled);
}

void TW_CALL TestConcaveConcave(void *)
{
    PTr.TestGetNodes(TVConcave,
                     TVConcave,
                     TraceDirect,
                     Emitter,Receiver,Disabled);
}

void TW_CALL TestConcaveFlat(void *)
{
    PTr.TestGetNodes(TVConcave,
                     TVFlat,
                     TraceDirect,
                     Emitter,Receiver,Disabled);
}

void TW_CALL TestFlatFlat(void *)
{
    PTr.TestGetNodes(TVFlat,
                     TVFlat,
                     TraceDirect,
                     Emitter,Receiver,Disabled);
}


void TW_CALL TestLoops(void *)
{
    PTr.TestGetNodes(TVInternal,
                     TVInternal,
                     TraceLoop,
                     Emitter,Receiver,Disabled);
}

//void TW_CALL SplitSupPatches(void *)
//{
//    PTr.SplitIntoSubPaths();
//    PTr.InitMetaMesh();
//    PTr.GetCurrCandidates(CurrCandidates);

//    ChosenCandidates.clear();
//    PTr.GetCurrChosen(ChosenCandidates);

//    ChosenIsLoop.clear();
//    PTr.GetCurrChosenIsLoop(ChosenIsLoop);

//    PTr.GetUnsatisfied(UnsatisfiedNodes);

//    PTr.GetVisualCornersPos(PatchCornerPos);
//}

//void TW_CALL SplitInIntervals(void *)
//{
//    PTr.SplitIntoIntervals();
//    UpdateVisualNodes();
////    PTr.InitMetaMesh();
////    PTr.GetCurrCandidates(CurrCandidates);

////    ChosenCandidates.clear();
////    PTr.GetCurrChosen(ChosenCandidates);

////    ChosenIsLoop.clear();
////    PTr.GetCurrChosenIsLoop(ChosenIsLoop);

////    PTr.GetUnsatisfied(UnsatisfiedNodes);

////    PTr.GetVisualCornersPos(PatchCornerPos);
//}

void TW_CALL SaveData(void *)
{
    SaveAllData(PTr,pathProject,CurrNum,subdivide_when_save,has_original_faces);
    SaveCSV(PTr,pathProject,CurrNum);
}

void  ProcessAllBatch()
{
    std::cout<<"TEST BATCH"<<std::endl;

    InitStructures();
    //std::cout<<"DE BOIA 4:"<<PTr.EDirTable.ConvexV.size()<<std::endl;

    //RecursiveProcess<TraceMesh>(PTr,Drift, add_only_needed,final_removal,);//,interleave_smooth);
    RecursiveProcess<TracerType>(PTr,Drift, add_only_needed,final_removal,true,meta_mesh_collapse,force_split);
    //RecursiveProcess<TraceMesh>(PTr,Drift,true,true,true);
    CurrPatchMode=CMPatchCol;
    //    PTr.BatchRemoval();
    //    PTr.FixValences();
    CurrPatchMode=CMPatchCol;
    drawField=false;
    drawSharpF=false;
    drawSing=false;
    UpdateVisualNodes();
    PTr.SmoothPatches();
    SaveAllData(PTr,pathProject,CurrNum,subdivide_when_save,has_original_faces);

    if (save_csv)
        SaveCSV(PTr,pathProject,CurrNum);

    if (save_setup)
        SaveSetupFile(pathProject,CurrNum);
}

void TW_CALL AllProcess(void *)
{
    ProcessAllBatch();
}

void InitLoopBar(QWidget *w)
{
    (void)w;
    bar = TwNewBar("LoopTool");
    TwDefine("LoopTool size='700 1500' ");
    TwDefine("LoopTool position='40 40' ");

    TwCopyCDStringToClientFunc (CopyCDStringToClient);


    TwAddButton(bar,"Reload All",ReloadAll,0," label='Reload All' ");

    TwEnumVal drawmodes[4] = { {GLW::DMSmooth, "Smooth"}, {GLW::DMPoints, "Per Points"},{GLW::DMFlatWire, "FlatWire"},{GLW::DMFlat, "Flat"}};
    // Create a type for the enum shapeEV
    TwType drawMode = TwDefineEnum("DrawMode", drawmodes, 4);
    TwAddVarRW(bar, "Draw Mode", drawMode, &drawmode, " keyIncr='<' keyDecr='>' help='Change draw mode.' ");

    TwEnumVal patchcolmodes[9] = { {CMPatchNone, "None"}, {CMPatchCol, "Per Patch"},
                                   {CMPatchValence, "Valence"},{CMPatchTopo, "Topology"},
                                   {CMPatchLengthDist,"Length Dist"},{CMPatchLenghtVar,"Length Var"},
                                   {CMArapDistortion,"Arap Dist"},{CMCClarkability,"Catmull Clarkability"},
                                   {CMExpValence,"Expected Valence"}
                                 };
    // Create a type for the enum shapeEV
    TwType patchColMode = TwDefineEnum("PatchColMode", patchcolmodes, 9);
    TwAddVarRW(bar, "Patch Col Mode", patchColMode, &CurrPatchMode, " keyIncr='<' keyDecr='>' help='Change col mode.' ");


    TwAddSeparator(bar,NULL,NULL);

    TwAddVarRW(bar,"DrawSharpF",TW_TYPE_BOOLCPP,&drawSharpF,"label='Draw Sharp Features'");
    TwAddVarRW(bar,"DrawSing",TW_TYPE_BOOLCPP,&drawSing,"label='Draw Singularities'");
    TwAddVarRW(bar,"DrawField",TW_TYPE_BOOLCPP,&drawField,"label='Draw Field'");

    TwAddSeparator(bar,NULL,NULL);
    TwAddVarRW(bar,"DrawTwins",TW_TYPE_BOOLCPP,&drawTwins,"label='Draw Twins'");

    TwAddVarRW(bar,"DrawNarrowEmit",TW_TYPE_BOOLCPP,&drawNarrowEmitters,"label='Draw Narrow Emitters'");
    TwAddVarRW(bar,"DrawNarrowReceive",TW_TYPE_BOOLCPP,&drawNarrowReceivers,"label='Draw Narrow Receivers'");


    TwAddVarRW(bar,"DrawPaths",TW_TYPE_BOOLCPP,&drawPaths,"label='Draw Paths'");
    TwAddVarRW(bar,"DrawPathNodes",TW_TYPE_BOOLCPP,&drawPathNodes,"label='Draw Path Nodes'");

    TwAddVarRW(bar,"DrawConcaveEmit",TW_TYPE_BOOLCPP,&drawConcaveEmitters,"label='Draw Concave Emitters'");
    TwAddVarRW(bar,"DrawConcaveReceive",TW_TYPE_BOOLCPP,&drawConcaveReceivers,"label='Draw Concave Receivers'");

    TwAddVarRW(bar,"DrawFlatEmit",TW_TYPE_BOOLCPP,&drawFlatEmitters,"label='Draw Flat Emitters'");
    TwAddVarRW(bar,"DrawFlatReceive",TW_TYPE_BOOLCPP,&drawFlatReceivers,"label='Draw Flat Receivers'");

    TwAddVarRW(bar,"DrawChoosenEmit",TW_TYPE_BOOLCPP,&drawChoosenEmitters,"label='Draw Choosen Emitters'");
    TwAddVarRW(bar,"DrawChoosenReceive",TW_TYPE_BOOLCPP,&drawChoosenReceivers,"label='Draw Choosen Receivers'");

    TwAddVarRW(bar,"DrawInvalidated",TW_TYPE_BOOLCPP,&drawInvalidated,"label='Draw Invalidated'");

    TwAddVarRW(bar,"DrawUnsatisfied",TW_TYPE_BOOLCPP,&drawUnsatisfied,"label='Draw Unsatisfied'");

    TwAddVarRW(bar,"DrawCorners",TW_TYPE_BOOLCPP,&drawCorners,"label='Draw Corners'");

    TwAddVarRW(bar,"DrawMetamesh",TW_TYPE_BOOLCPP,&drawMetaMesh,"label='Draw MetaMesh'");

    TwAddVarRW(bar,"DrawProblematics",TW_TYPE_BOOLCPP,&drawProblematicsOnly,"label='Draw Problematic'");

    TwAddVarRW(bar,"VIsTraces",TW_TYPE_INT32,&VisTraces,"label='Visualize Until'");


    TwAddSeparator(bar,NULL,NULL);

    TwAddVarRW(bar,"Drift",TW_TYPE_DOUBLE,&Drift,"label='Drift'");
    TwAddButton(bar,"InitGraph",InitGraph,0," label='Init Graph' ");
    TwAddButton(bar,"JoinNarrow",JoinNarrow,0," label='Trace Narrow' ");
    TwAddButton(bar,"JoinConcave",JoinConcave,0," label='Trace Concave' ");
    TwAddButton(bar,"TraceLoops",AddLoops,0," label='Trace Loops' ");
    TwAddButton(bar,"TraceBorders",TraceBorder,0," label='Trace Borders' ");
    TwAddButton(bar,"BatchProcess",BatchProcess,0," label='Batch Process' ");
    //TwAddButton(bar,"Split",SplitSupPatches,0," label='Split sub' ");
    TwAddButton(bar,"BatchRemoval",BatchRemoval,0," label='Batch Removal' ");
    TwAddButton(bar,"SmoothPaths",SmoothPathes,0," label='Smooth Paths' ");
    //TwAddButton(bar,"SmoothParam",SmoothParam,0," label='Smooth Param' ");


    TwAddSeparator(bar,NULL,NULL);

    TwAddVarRW(bar,"Sample Ratio",TW_TYPE_DOUBLE,&PTr.sample_ratio,"label='Sample Ratio'");



    TwAddVarRW(bar,"SplitOnRemove",TW_TYPE_BOOLCPP,&PTr.split_on_removal,"label='Split on Remove'");
    TwAddVarRW(bar,"MaxValence",TW_TYPE_INT32,&PTr.MaxVal,"label='Max Valence'");
    TwAddVarRW(bar,"CCLarkability",TW_TYPE_DOUBLE,&PTr.CClarkability,"label='CCLarkability'");
    TwAddVarRW(bar,"ConcaveNeed",TW_TYPE_INT32,&PTr.Concave_Need,"label='Concave Need'");
    //TwAddVarRW(bar,"MaxValence",TW_TYPE_INT32,&PTr.MaxVal,"label='Max Valence'");
    TwAddVarRW(bar,"MatchVal",TW_TYPE_BOOLCPP,&PTr.match_valence,"label='Match Valence'");
    TwAddVarRW(bar,"AddNeed",TW_TYPE_BOOLCPP,&add_only_needed,"label='Add Only need'");
    TwAddVarRW(bar,"MetaMeshColl",TW_TYPE_BOOLCPP,&meta_mesh_collapse,"label='Meta Mesh Collapse'");

    //TwAddVarRW(bar,"InterRem",TW_TYPE_BOOLCPP,&interleave_removal,"label='Interleave removal'");
    //TwAddVarRW(bar,"InterSmth",TW_TYPE_BOOLCPP,&interleave_smooth,"label='Interleave smooth'");
    TwAddVarRW(bar,"FinalRem",TW_TYPE_BOOLCPP,&final_removal,"label='Final removal'");
    TwAddVarRW(bar,"ForceSplit",TW_TYPE_BOOLCPP,&force_split,"label='Force split'");
    TwAddButton(bar,"RecursiveProcess",RecursiveProcess,0," label='Recursive Process' ");
    TwAddButton(bar,"Subdivide",SubdividePatches,0," label='Subdivide Patches' ");

  //    TwAddSeparator(bar,NULL,NULL);

    TwAddButton(bar,"Subdivide",SubdividePatches,0," label='Subdivide Patches' ");
//    TwAddButton(bar,"ParametrizePatches",ParametrizePatches,0," label='Parametrize Patches' ");


    //TwAddSeparator(bar,NULL,NULL);

    //TwAddButton(bar,"TestInit",TestInit,0," label='TestInit' ");
    //TwAddButton(bar,"Split",SplitSupPatches,0," label='Split sub' ");
    //TwAddButton(bar,"Split Int",SplitInIntervals,0," label='Split in intervals' ");
    //TwAddButton(bar,"OneStepRem",RemoveIteration,0," label='Remove Ite step' ");
    TwAddButton(bar,"ParametrizePatches",ParametrizePatches,0," label='Parametrize Patches' ");


    TwAddSeparator(bar,NULL,NULL);

    TwAddVarRW(bar,"SubdSave",TW_TYPE_BOOLCPP,&subdivide_when_save,"label='Subdivide When Save'");
    TwAddButton(bar,"SaveData",SaveData,0," label='Save Data' ");


    TwAddSeparator(bar,NULL,NULL);

    TwAddButton(bar,"AllProcess",AllProcess,0," label='All Process' ");

}

void InitBar(QWidget *w)
{
    InitLoopBar(w);
}


GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
    filename=0;
    //hasToPick=false;
    InitBar(this);
    //LoadSetup();
    LoadSetupFile(std::string("basic_setup.txt"));
    FindCurrentNum();
    LoadAll();
    if (batch_process)
    {
        ProcessAllBatch();
        exit(0);
    }
}

void GLWidget::initializeGL ()
{
    //    FieldParam.alpha_curv=0.2;
    //    FieldParam.curvRing=4;
    glewInit();
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
    initializeGL();
}

void GLWidget::paintGL ()
{

    if (CurrPatchMode!=OldPatchMode)
    {
        UpdatePatchColor();
        OldPatchMode=CurrPatchMode;
    }

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
    //glPushMatrix();
    glWrap.m = &mesh;
    glWrapProblem.m=&problematic_mesh;
    if(mesh.vert.size()>0)
    {
        vcg::glScale(2.0f/mesh.bbox.Diag());
        glTranslate(-mesh.bbox.Center());

        vcg::glColor(vcg::Color4b(220,220,220,255));
        //glWrap.Draw(GLW::DrawMode(drawmode),GLW::CMNone,GLW::TMNone);

        if (!drawProblematicsOnly)
            glWrap.Draw(GLW::DrawMode(drawmode),GLW::CMPerFace,GLW::TMNone);
        else
        {
            glWrap.Draw(GLW::DMWire,GLW::CMNone,GLW::TMNone);
            glWrapProblem.Draw(GLW::DrawMode(drawmode),GLW::CMPerFace,GLW::TMNone);
        }

        if (drawField)
        {
            if (!drawProblematicsOnly)
                vcg::GLField<TraceMesh>::GLDrawFaceField(mesh,false,false,0.005);
            else
                vcg::GLField<TraceMesh>::GLDrawFaceField(problematic_mesh,false,false,0.005);
        }

        if (drawSing)
        {
            if (!drawProblematicsOnly)
                vcg::GLField<TraceMesh>::GLDrawSingularity(mesh);
            else
                vcg::GLField<TraceMesh>::GLDrawFaceField(problematic_mesh,false,false,0.005);
        }

        if (drawTwins)
            GLGraph.GLDrawTwinsConnections();
        if (drawConcaveEmitters)
            GLGraph.GLDrawNodes(ConcaveEmittersNode,mesh.bbox.Diag()*0.005);
        if (drawConcaveReceivers)
            GLGraph.GLDrawNodes(ConcaveReceiversNode,mesh.bbox.Diag()*0.005);
        if (drawFlatEmitters)
            GLGraph.GLDrawNodes(FlatEmittersNode,mesh.bbox.Diag()*0.005,false,10,true);
        if (drawChoosenEmitters)
            GLGraph.GLDrawNodes(ChoosenEmittersNode,mesh.bbox.Diag()*0.005);
        if (drawFlatReceivers)
            GLGraph.GLDrawNodes(FlatReceiversNode,mesh.bbox.Diag()*0.005);
        if (drawNarrowEmitters)
            GLGraph.GLDrawNodes(NarrowEmittersNode,mesh.bbox.Diag()*0.005);
        if (drawNarrowReceivers)
            GLGraph.GLDrawNodes(NarrowReceiversNode,mesh.bbox.Diag()*0.005);
        if (drawChoosenReceivers)
            GLGraph.GLDrawNodes(ChoosenReceiversNode,mesh.bbox.Diag()*0.005);

        if (drawInvalidated)
            GLGraph.GLDrawNonActiveNodes(mesh.bbox.Diag()*0.001);
        if (drawSharpF)
            MeshDrawing<TraceMesh>::GLDrawSharpEdges(mesh);
            //mesh.GLDrawSharpEdges();
        if (drawUnsatisfied)
            GLGraph.GLDrawNodes(UnsatisfiedNodes,mesh.bbox.Diag()*0.001,true,30);

        if (testDrawEmitter)
            GLGraph.GLDrawNodes(Emitter,mesh.bbox.Diag()*0.002);
        if (testdrawReceiver)
            GLGraph.GLDrawNodes(Receiver,mesh.bbox.Diag()*0.002);
        if (testdrawDisables)
            GLGraph.GLDrawNodes(Disabled,mesh.bbox.Diag()*0.002);


        if (drawCorners)
            GLGraph.GLDrawPoints(PatchCornerPos,10,vcg::Color4b(255,0,255,255));

        if (drawPaths)
            GLGraph.GLDrawPaths(ChosenCandidates,ChosenIsLoop,mesh.bbox.Diag()*0.01,drawPathNodes,VisTraces);//,true);

        //GLGraph.GLDrawPaths(DiscardedCandidates,DiscardedIsLoop,mesh.bbox.Diag()*0.01,true);
        //GLGraph.GLDrawNodes(TraceableFlatNode,mesh.bbox.Diag()*0.002);

        if (drawMetaMesh)
            GLMetaMesh<TraceMesh>::GLDraw(PTr.MMesh);
            //PTr.GLDraweMetaMesh();
        //GLGraph.GLDrawSingNodes(mesh.bbox.Diag()*0.002);
    }

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
    (void)e;
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
