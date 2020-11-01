/****************************************************************************
**
** Copyright (C) 2011 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
**     the names of its contributors may be used to endorse or promote
**     products derived from this software without specific prior written
**     permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
** $QT_END_LICENSE$
**
****************************************************************************/

#include "glwidget.h"
#include <wrap/qt/trackball.h>
#include <wrap/gl/picking.h>
#include <wrap/qt/anttweakbarMapper.h>
#include <wrap/gl/gl_field.h>
#include <QDir>
#include "tracing/GL_vert_field_graph.h"
#include "tracing/patch_tracer.h"
#include "tracing/tracer_interface.h"

TwBar *bar;
char * filename;/// filename of the mesh to load
CMesh mesh;     /// the active mesh instance
vcg::GlTrimesh<CMesh> glWrap;    /// the active mesh opengl wrapper
vcg::Trackball track;     /// the active manipulator
GLW::DrawMode drawmode=GLW::DMFlatWire;     /// the current drawmode

std::string pathM,pathF,pathS,pathProject;

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

bool drawNarrowCandidates=false;

//bool splitted=false;

bool batch_process=false;
bool has_features=false;

bool add_only_needed=true;
bool interleave_removal=true;
bool interleave_smooth=false;
bool final_removal=true;

int CurrNum=0;

std::vector<size_t> ConcaveEmittersNode,ConcaveReceiversNode,
FlatEmittersNode,FlatReceiversNode,
NarrowEmittersNode,NarrowReceiversNode,
ProblematicNodes,UnsatisfiedNodes,
ChoosenEmittersNode,ChoosenReceiversNode;

VertexFieldGraph<CMesh> VGraph(mesh);
GLVertGraph<CMesh> GLGraph(VGraph);

PatchTracer<CMesh> PTr(VGraph);

std::vector<std::vector<size_t> > CurrCandidates;
std::vector<bool> ChosenIsLoop;
std::vector<std::vector<size_t> > ChosenCandidates;
std::vector<typename CMesh::CoordType> PatchCornerPos;
CMesh::ScalarType Drift=100;

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
//    Drift=(CMesh::ScalarType)Driftf;
//    std::cout<<"DRIFT "<<Drift<<std::endl;
    fprintf(f,"Drift %f\n",Drift);

//    float SRatef;
//    fscanf(f,"Srate %f\n",&SRatef);
//    PTr.sample_ratio=(CMesh::ScalarType)SRatef;
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
//    PTr.CClarkability=(CMesh::ScalarType)CCbility;
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
    if (interleave_removal)
        fprintf(f,"InterleaveRem 1\n");
    else
        fprintf(f,"InterleaveRem 0\n");

//    fscanf(f,"InterleaveSmooth %d\n",&IntVar);
//    std::cout<<"INTERLEAVE SMOOTH "<<IntVar<<std::endl;
//    if (IntVar==0)
//        interleave_smooth=false;
//    else
//        interleave_smooth=true;
    if (interleave_smooth)
        fprintf(f,"InterleaveSmooth 1\n");
    else
        fprintf(f,"InterleaveSmooth 0\n");

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
    Drift=(CMesh::ScalarType)Driftf;
    std::cout<<"DRIFT "<<Drift<<std::endl;

    float SRatef;
    fscanf(f,"Srate %f\n",&SRatef);
    PTr.sample_ratio=(CMesh::ScalarType)SRatef;
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
    PTr.CClarkability=(CMesh::ScalarType)CCbility;
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

    fscanf(f,"InterleaveRem %d\n",&IntVar);
    std::cout<<"INTERLEAVE REMOVE "<<IntVar<<std::endl;
    if (IntVar==0)
        interleave_removal=false;
    else
        interleave_removal=true;

    fscanf(f,"InterleaveSmooth %d\n",&IntVar);
    std::cout<<"INTERLEAVE SMOOTH "<<IntVar<<std::endl;
    if (IntVar==0)
        interleave_smooth=false;
    else
        interleave_smooth=true;

    fscanf(f,"FinalRem %d\n",&IntVar);
    std::cout<<"FINAL REMOVE "<<IntVar<<std::endl;
    if (IntVar==0)
        final_removal=false;
    else
        final_removal=true;

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
    std::cout<<"basepath "<<BasePath.c_str()<<std::endl;
    QDir directory(BasePath.c_str());

    QFile f(pathM.c_str());
    QFileInfo fileInfo(f.fileName());
    QString filename(fileInfo.fileName());
    std::string NameFile=filename.toStdString();
    std::cout<<"namefile "<<NameFile.c_str()<<std::endl;
    NameFile.resize(NameFile.find_last_of("."));
    std::string Filter=NameFile+"_p*.obj";
    std::cout<<"filter "<<Filter.c_str()<<std::endl;

//    TestPath+="*.obj";
   QStringList projectFiles = directory.entryList(QStringList() <<Filter.c_str(),QDir::Files);
    CurrNum=0;

    foreach(QString filename, projectFiles) {
        int Num;
        std::string TestScan=NameFile+"_p%d.obj";
        sscanf (filename.toStdString().c_str(),TestScan.c_str(),&Num);
        CurrNum=std::max(CurrNum,(Num+1));
        std::cout<<"test "<<Num<<std::endl;
//        std::cout<<"test "<<filename.toStdString().c_str()<<std::endl;
    //do whatever you need to do
    }
}

//void LoadSetupFile(std::string path)
//{
//    FILE *f=fopen(path.c_str(),"rt");
//    assert(f!=NULL);
//    float Driftf;
//    fscanf(f,"Drift %f\n",&Driftf);
//    Drift=(CMesh::ScalarType)Driftf;
//    std::cout<<"DRIFT "<<Drift<<std::endl;

//    int IntVar=0;
//    fscanf(f,"Split %d\n",&IntVar);
//    std::cout<<"SPLIT "<<IntVar<<std::endl;
//    if (IntVar==0)
//        PTr.split_on_removal=false;
//    else
//        PTr.split_on_removal=true;

//    fscanf(f,"IncreaseValRem %d\n",&IntVar);
//    std::cout<<"INCREASE VAL "<<IntVar<<std::endl;
//    if (IntVar==0)
//        PTr.avoid_increase_valence=false;
//    else
//        PTr.avoid_increase_valence=true;

//    fscanf(f,"IrregularRem %d\n",&IntVar);
//    std::cout<<"IRR VAL "<<IntVar<<std::endl;
//    if (IntVar==0)
//        PTr.avoid_collapse_irregular=false;
//    else
//        PTr.avoid_collapse_irregular=true;

//    float MaxDistortionf;
//    fscanf(f,"DistortionL %f\n",&MaxDistortionf);
//    PTr.max_lenght_distortion=(CMesh::ScalarType)MaxDistortionf;
//    fclose(f);
//}

void UpdatePatchColor()
{
 if (CurrPatchMode==CMPatchNone)
     vcg::tri::UpdateColor<CMesh>::PerFaceConstant(mesh);
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

    ChosenIsLoop.clear();


    PTr.GetCurrChosenIsLoop(ChosenIsLoop);

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
//    for (size_t i=0;i<UnsatisfiedNodes.size();i++)
//        std::cout<<"Unsatisfied "<<UnsatisfiedNodes[i]<<std::endl;
    PTr.GetVisualCornersPos(PatchCornerPos);
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


    PreProcessMesh(mesh);

    VGraph.Init();//SingPos);


    GLGraph.InitDisplacedPos();

    PTr.Init(Drift);
    //std::cout<<"Here"<<std::endl;
    UpdateVisualNodes();
}

void TW_CALL InitGraph(void *)
{
    InitStructures();

    drawField=false;
    drawSharpF=false;
    drawSing=false;
}

void TW_CALL JoinNarrow(void *)
{
    PTr.JoinNarrow();
    UpdateVisualNodes();
}

void TW_CALL JoinConcave(void *)
{
    PTr.JoinConcave();
    UpdateVisualNodes();
}

void TW_CALL AddLoops(void *)
{
    PTr.TraceLoops();
    UpdateVisualNodes();
}

void TW_CALL TraceBorder(void *)
{
    PTr.JoinBoundaries();
    UpdateVisualNodes();
}

void TW_CALL SmoothPathes(void *)
{
    PTr.SmoothPatches();
    PTr.GetVisualCornersPos(PatchCornerPos);
}


//void TW_CALL RemovePath(void *)
//{
//    PTr.RemovePaths();
//    CurrCandidates.clear();
//    PTr.GetCurrCandidates(CurrCandidates);

//    ChosenCandidates.clear();
//    PTr.GetCurrChosen(ChosenCandidates);

//    ChosenIsLoop.clear();
//    PTr.GetCurrChosenIsLoop(ChosenIsLoop);

//    PTr.GetUnsatisfied(UnsatisfiedNodes);

//    PTr.GetVisualCornersPos(PatchCornerPos);

////    PTr.GetChoosenEmitters(ChoosenEmittersNode);
////    PTr.GetChoosenReceivers(ChoosenReceiversNode);
//}


void TW_CALL BatchProcess(void *)
{
    InitStructures();
    //PTr.BatchProcess();
    PTr.BatchAddLoops(true,false,false,true,false);
    CurrPatchMode=CMPatchCol;
//    CurrCandidates.clear();
//    PTr.GetCurrCandidates(CurrCandidates);

//    ChosenCandidates.clear();
//    PTr.GetCurrChosen(ChosenCandidates);
//    ChosenIsLoop.clear();
//    PTr.GetCurrChosenIsLoop(ChosenIsLoop);

//    PTr.GetUnsatisfied(UnsatisfiedNodes);

    //PTr.GetCornersPos(PatchCornerPos);

    drawField=false;
    drawSharpF=false;
    drawSing=false;

    UpdateVisualNodes();
//    PTr.GetConcaveEmitters(ConcaveEmittersNode);
//    PTr.GetConcaveReceivers(ConcaveReceiversNode);
//    PTr.GetFlatEmitters(FlatEmittersNode);
//    PTr.GetFlatReceivers(FlatReceiversNode);
//    PTr.GetNarrowActiveEmitters(NarrowEmittersNode);
//    PTr.GetNarrowActiveReceivers(NarrowReceiversNode);
//    PTr.GetChoosenEmitters(ChoosenEmittersNode);
//    PTr.GetChoosenReceivers(ChoosenReceiversNode);
//    PTr.GetUnsatisfied(UnsatisfiedNodes);
//    PTr.GetCornersPos(PatchCornerPos);
}

//void TW_CALL IterativeBatch(void *)
//{
//    InitStructures();
//    //RecursiveProcess<CMesh>(PTr,Drift);
//    RecursiveProcess3<CMesh>(PTr,Drift);
////    CurrPatchMode=CMPatchCol;
////    drawField=false;
////    drawSharpF=false;
////    drawSing=false;
////    UpdateVisualNodes();
//}

void TW_CALL RecursiveProcess(void *)
{
    InitStructures();
    //RecursiveProcess<CMesh>(PTr,Drift);
    RecursiveProcess<CMesh>(PTr,Drift, add_only_needed,interleave_removal,final_removal,interleave_smooth);
    CurrPatchMode=CMPatchCol;
    drawField=false;
    drawSharpF=false;
    drawSing=false;
    UpdateVisualNodes();
//    ChosenCandidates.clear();
//    PTr.GetCurrChosen(ChosenCandidates);
//    ChosenIsLoop.clear();
//    PTr.GetCurrChosenIsLoop(ChosenIsLoop);
//    PTr.GetUnsatisfied(UnsatisfiedNodes);
}

void TW_CALL ParametrizePatches(void *)
{
    PTr.ComputePatchesUV();
//    InitStructures();
//    //RecursiveProcess<CMesh>(PTr,Drift);
//    RecursiveProcess<CMesh>(PTr,Drift, add_only_needed,interleave_removal,final_removal);
//    CurrPatchMode=CMPatchCol;
//    drawField=false;
//    drawSharpF=false;
//    drawSing=false;
//    UpdateVisualNodes();
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
    PTr.BatchRemoval(false);
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

    if (!has_features)return;

    bool loadedFeatures=mesh.LoadSharpFeatures(pathS);
    if (!loadedFeatures){
        std::cout<<"*** ERROR LOADING FEATURES ***"<<std::endl;
        exit(0);
    }
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

void TW_CALL SplitSupPatches(void *)
{
    PTr.SplitIntoSubPaths();
    PTr.GetCurrCandidates(CurrCandidates);

    ChosenCandidates.clear();
    PTr.GetCurrChosen(ChosenCandidates);

    ChosenIsLoop.clear();
    PTr.GetCurrChosenIsLoop(ChosenIsLoop);

    PTr.GetUnsatisfied(UnsatisfiedNodes);

    PTr.GetVisualCornersPos(PatchCornerPos);
}


void TW_CALL SaveData(void *)
{
    SaveAllData(PTr,pathProject,CurrNum);
    SaveCSV(PTr,pathProject,CurrNum);
}

void  ProcessAllBatch()
{
    InitStructures();
    RecursiveProcess<CMesh>(PTr,Drift, add_only_needed,interleave_removal,final_removal,interleave_smooth);
    //RecursiveProcess<CMesh>(PTr,Drift,true,true,true);
    CurrPatchMode=CMPatchCol;
    PTr.BatchRemoval();
    PTr.FixValences();
    CurrPatchMode=CMPatchCol;
    drawField=false;
    drawSharpF=false;
    drawSing=false;
    UpdateVisualNodes();
    PTr.SmoothPatches();
    SaveAllData(PTr,pathProject,CurrNum);
    SaveCSV(PTr,pathProject,CurrNum);
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

    TwAddSeparator(bar,NULL,NULL);

    TwAddVarRW(bar,"Drift",TW_TYPE_DOUBLE,&Drift,"label='Drift'");
    TwAddButton(bar,"InitGraph",InitGraph,0," label='Init Graph' ");
    TwAddButton(bar,"JoinNarrow",JoinNarrow,0," label='Trace Narrow' ");
    TwAddButton(bar,"JoinConcave",JoinConcave,0," label='Trace Concave' ");
    TwAddButton(bar,"TraceLoops",AddLoops,0," label='Trace Loops' ");
    TwAddButton(bar,"TraceBorders",TraceBorder,0," label='Trace Borders' ");
    TwAddButton(bar,"BatchProcess",BatchProcess,0," label='Batch Process' ");
    TwAddButton(bar,"Split",SplitSupPatches,0," label='Split sub' ");
    TwAddButton(bar,"BatchRemoval",BatchRemoval,0," label='Batch Removal' ");
    TwAddButton(bar,"SmoothPaths",SmoothPathes,0," label='Smooth Paths' ");

    TwAddSeparator(bar,NULL,NULL);

    TwAddVarRW(bar,"Sample Ratio",TW_TYPE_DOUBLE,&PTr.sample_ratio,"label='Sample Ratio'");



    TwAddVarRW(bar,"SplitOnRemove",TW_TYPE_BOOLCPP,&PTr.split_on_removal,"label='Split on Remove'");
    TwAddVarRW(bar,"MaxValence",TW_TYPE_INT32,&PTr.MaxVal,"label='Max Valence'");
    TwAddVarRW(bar,"CCLarkability",TW_TYPE_DOUBLE,&PTr.CClarkability,"label='CCLarkability'");
    TwAddVarRW(bar,"ConcaveNeed",TW_TYPE_INT32,&PTr.Concave_Need,"label='Concave Need'");
    TwAddVarRW(bar,"MaxValence",TW_TYPE_INT32,&PTr.MaxVal,"label='Max Valence'");
    TwAddVarRW(bar,"AddNeed",TW_TYPE_BOOLCPP,&add_only_needed,"label='Add Only need'");
    TwAddVarRW(bar,"InterRem",TW_TYPE_BOOLCPP,&interleave_removal,"label='Interleave removal'");
    TwAddVarRW(bar,"InterSmth",TW_TYPE_BOOLCPP,&interleave_smooth,"label='Interleave smooth'");
    TwAddVarRW(bar,"FinalRem",TW_TYPE_BOOLCPP,&final_removal,"label='Final removal'");
    TwAddButton(bar,"RecursiveProcess",RecursiveProcess,0," label='Recursive Process' ");
    TwAddButton(bar,"Subdivide",SubdividePatches,0," label='Subdivide Patches' ");
    TwAddButton(bar,"ParametrizePatches",ParametrizePatches,0," label='Parametrize Patches' ");

    TwAddButton(bar,"SaveData",SaveData,0," label='Save Data' ");


    TwAddSeparator(bar,NULL,NULL);

    TwAddButton(bar,"AllProcess",AllProcess,0," label='All Process' ");



//    TwAddSeparator(bar,NULL,NULL);
//    TwAddButton(bar,"TestNarrowNarrow",TestNarrowNarrow,0," label='TestNarrowNarrow' ");
//    TwAddButton(bar,"TestNarrowConcave",TestNarrowConcave,0," label='TestNarrowConcave' ");
//    TwAddButton(bar,"TestNarrowFlat",TestNarrowFlat,0," label='TestNarrowFlat' ");
////    TwAddButton(bar,"TestNarrowChosen",TestNarrowChosen,0," label='TestNarrowChosen' ");

//    TwAddButton(bar,"TestConcaveConcave",TestConcaveConcave,0," label='TestConcaveConcave' ");
//    TwAddButton(bar,"TestConcaveFlat",TestConcaveFlat,0," label='TestConcaveFlat' ");
////    TwAddButton(bar,"TestConcaveChosen",TestConcaveChosen,0," label='TestConcaveChosen' ");

//    TwAddButton(bar,"TestFlatFlat",TestFlatFlat,0," label='TestFlatFlat' ");
////    TwAddButton(bar,"TestFlatChosen",TestFlatChosen,0," label='TestFlatChosen' ");

////    TwAddButton(bar,"TestChosenChosen",TestChosenChosen,0," label='TestChosenChosen' ");
//    TwAddButton(bar,"TestLoops",TestLoops,0," label='TestLoops' ");

//    TwAddVarRW(bar,"testDrawEmitter",TW_TYPE_BOOLCPP,&testDrawEmitter,"label='testDrawEmitter'");
//    TwAddVarRW(bar,"testdrawReceiver",TW_TYPE_BOOLCPP,&testdrawReceiver,"label='testdrawReceiver'");
//    TwAddVarRW(bar,"testdrawDisables",TW_TYPE_BOOLCPP,&testdrawDisables,"label='testdrawDisables'");

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
    if(mesh.vert.size()>0)
    {
        vcg::glScale(2.0f/mesh.bbox.Diag());
        glTranslate(-mesh.bbox.Center());

        vcg::glColor(vcg::Color4b(220,220,220,255));
        //glWrap.Draw(GLW::DrawMode(drawmode),GLW::CMNone,GLW::TMNone);

        glWrap.Draw(GLW::DrawMode(drawmode),GLW::CMPerFace,GLW::TMNone);

        if (drawField)//&&(has_field))
            vcg::GLField<CMesh>::GLDrawFaceField(mesh,false,false,0.007);

        if (drawSing)//&&(has_field))
            vcg::GLField<CMesh>::GLDrawSingularity(mesh);

        if (drawTwins)
            GLGraph.GLDrawTwinsConnections();
        if (drawConcaveEmitters)
            GLGraph.GLDrawNodes(ConcaveEmittersNode,mesh.bbox.Diag()*0.005);
        if (drawConcaveReceivers)
            GLGraph.GLDrawNodes(ConcaveReceiversNode,mesh.bbox.Diag()*0.005);
        if (drawFlatEmitters)
            GLGraph.GLDrawNodes(FlatEmittersNode,mesh.bbox.Diag()*0.005);
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
            mesh.GLDrawSharpEdges();
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
            GLGraph.GLDrawPaths(ChosenCandidates,ChosenIsLoop,mesh.bbox.Diag()*0.01,drawPathNodes);

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
