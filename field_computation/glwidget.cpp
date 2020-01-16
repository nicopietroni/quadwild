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

#include <GL/glew.h>
#include <QMouseEvent>

#include <math.h>
#include "glwidget.h"

#include <wrap/igl/smooth_field.h>
#include <triangle_mesh_type.h>
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

std::string pathM="";

vcg::Trackball track;//the active manipulator

bool drawfield=false;

MyTriMesh tri_mesh;

//MyTriMesh remeshed_mesh;
//MyTriMesh original_mesh;

TwBar *barQuad;

vcg::GlTrimesh<MyTriMesh> glWrap;
vcg::GLField<MyTriMesh> glField;

typedef typename MyTriMesh::ScalarType ScalarType;
typedef typename MyTriMesh::CoordType CoordType;

int Iterations;
ScalarType EdgeStep;
ScalarType Multiplier=2;
ScalarType SharpFactor=6;
//ScalarType sharp_feature_thr=45;

int xMouse,yMouse;

vcg::GridStaticPtr<MyTriMesh::FaceType,MyTriMesh::ScalarType> Gr;

typedef FieldSmoother<MyTriMesh> FieldSmootherType;
FieldSmootherType::SmoothParam FieldParam;

AutoRemesher<MyTriMesh>::Params RemPar;
bool do_batch=false;
//size_t ErodeDilateSteps=4;

int remesher_iterations=15;
ScalarType remesher_aspect_ratio=0.3;
int remesher_termination_delta=10000;

ScalarType sharp_feature_thr=35;
int feature_erode_dilate=4;

void TW_CALL AutoRemesh(void *)
{

   tri_mesh.UpdateDataStructures();
//   vcg::tri::IsotropicRemeshing<MyTriMesh>::Params par;
//   par.minLength=tri_mesh.bbox.Diag()*0.005;
//   par.maxLength=tri_mesh.bbox.Diag()*0.01;
//   vcg::tri::IsotropicRemeshing<MyTriMesh>::Do(tri_mesh,par);

   RemPar.iterations   = remesher_iterations;
   RemPar.targetAspect = remesher_aspect_ratio;
   RemPar.targetDeltaFN= remesher_termination_delta;
   RemPar.userSelectedCreases = true;
   RemPar.surfDistCheck = true;

   std::shared_ptr<MyTriMesh> clean = AutoRemesher<MyTriMesh>::CleanMesh(tri_mesh,true);

   std::shared_ptr<MyTriMesh> ret=AutoRemesher<MyTriMesh>::Remesh(*clean,RemPar);
   tri_mesh.Clear();
   vcg::tri::Append<MyTriMesh,MyTriMesh>::Mesh(tri_mesh,(*ret));
   vcg::tri::Clean<MyTriMesh>::RemoveUnreferencedVertex(tri_mesh);
   vcg::tri::Allocator<MyTriMesh>::CompactEveryVector(tri_mesh);
   tri_mesh.UpdateDataStructures();
}

void TW_CALL InitSharpFeatures(void *)
{
    tri_mesh.UpdateDataStructures();
	tri_mesh.InitSharpFeatures(sharp_feature_thr);

}

void TW_CALL RefineIfNeeded(void *)
{
    tri_mesh.RefineIfNeeded();
    Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
}

void TW_CALL SmoothField(void *)
{
    tri_mesh.SmoothField(FieldParam);
    drawfield=true;

    //Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
}

void SaveAllData()
{
    std::string projM=pathM;
    size_t indexExt=projM.find_last_of(".");
    projM=projM.substr(0,indexExt);
    std::string meshName=projM+std::string("_rem.obj");
    std::string fieldName=projM+std::string("_rem.rosy");
    std::string sharpName=projM+std::string("_rem.sharp");
    std::cout<<"Saving Mesh TO:"<<meshName.c_str()<<std::endl;
    std::cout<<"Saving Field TO:"<<fieldName.c_str()<<std::endl;
    std::cout<<"Saving Sharp TO:"<<sharpName.c_str()<<std::endl;
    tri_mesh.SaveTriMesh(meshName.c_str());
    tri_mesh.SaveField(fieldName.c_str());
    tri_mesh.SaveSharpFeatures(sharpName.c_str());
    //tri_mesh.SaveTriMesh()
}

void TW_CALL SaveData(void *)
{
    SaveAllData();
}

void TW_CALL ErodeDilateFeatureStep(void *)
{
	tri_mesh.ErodeDilate(feature_erode_dilate);
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

    TwAddButton(barQuad,"SetSharp",InitSharpFeatures,0,"label='InitSharp'");

	TwAddVarRW(barQuad,"ErodeDilSteps",TW_TYPE_INT32,&feature_erode_dilate,"label='ErodeDilateSteps'");

    TwAddButton(barQuad,"Erode Dilate",ErodeDilateFeatureStep,0,"label='ErodeDilateSharp'");

    TwAddButton(barQuad,"AutoRemesh",AutoRemesh,0,"label='AutoRemesh'");

    TwAddButton(barQuad,"Refine",RefineIfNeeded,0,"label='Refine if needed'");

    TwAddVarRW(barQuad,"Alpha",TW_TYPE_DOUBLE, &FieldParam.alpha_curv," label='Alpha Curvature'");
    TwAddVarRW(barQuad,"HardCT",TW_TYPE_DOUBLE, &FieldParam.curv_thr," label='Hard Curv Thr'");
    TwAddVarRW(barQuad,"CurvRing",TW_TYPE_INT32,&FieldParam.curvRing,"label='Curvature Ring'");

    TwEnumVal smoothmodes[2] = {
        {vcg::tri::SMMiq,"MIQ"},
        {vcg::tri::SMNPoly,"NPoly"}
    };
    TwType smoothMode = TwDefineEnum("SmoothMode", smoothmodes, 2);
    TwAddVarRW(barQuad, "Smooth Mode", smoothMode, &FieldParam.SmoothM," label='Smooth Mode' ");

    TwAddButton(barQuad,"ComputeField",SmoothField,0,"label='Compute Field'");

    TwAddButton(barQuad,"SaveData",SaveData,0,"label='Save Data'");

}

void BatchProcess ()
{
        //SHARP FEATURE
	RemPar.iterations   = remesher_iterations;
	RemPar.targetAspect = remesher_aspect_ratio;
	RemPar.targetDeltaFN= remesher_termination_delta;
	RemPar.userSelectedCreases = true;
	RemPar.surfDistCheck = true;

	std::cout << "[fieldComputation] Mesh cleaning..." << std::endl;
	std::shared_ptr<MyTriMesh> clean = AutoRemesher<MyTriMesh>::CleanMesh(tri_mesh, false);

	std::cout << "[fieldComputation] Feature Extraction..." << std::endl;
	clean->UpdateDataStructures();
	clean->InitSharpFeatures(sharp_feature_thr);
	clean->ErodeDilate(feature_erode_dilate);

	//REMESH
	std::cout << "[fieldComputation] Initial Remeshing..." << std::endl;
	std::shared_ptr<MyTriMesh> ret=AutoRemesher<MyTriMesh>::Remesh(*clean,RemPar);
	{	
	    	std::string projM=pathM;
    		size_t indexExt=projM.find_last_of(".");
    		projM=projM.substr(0,indexExt);
		std::string meshName=projM+std::string("_remeshed.obj");
		ret->SaveTriMesh(meshName.c_str());
	}

	AutoRemesher<MyTriMesh>::SplitNonManifold(*ret);

    tri_mesh.Clear();
    vcg::tri::Append<MyTriMesh,MyTriMesh>::Mesh(tri_mesh,(*ret));
    tri_mesh.UpdateDataStructures();
    
    //REFINE IF NEEDED
    tri_mesh.InitSharpFeatures(sharp_feature_thr);
    tri_mesh.ErodeDilate(feature_erode_dilate);
    tri_mesh.RefineIfNeeded();

    //FIELD SMOOTH
    bool UseNPoly=tri_mesh.SufficientFeatures(SharpFactor);
    if (UseNPoly)
    {
        std::cout<<"Using NPoly"<<std::endl;
        FieldParam.SmoothM=SMNPoly;
    }
    else
    {
        std::cout<<"Using Comiso"<<std::endl;
        FieldParam.SmoothM=SMMiq;
        FieldParam.alpha_curv=0.3;
    }

    std::cout << "[fieldComputation] Smooth Field Computation..." << std::endl;
    tri_mesh.SmoothField(FieldParam);

    //SAVE
    SaveAllData();
}

GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
    hasToPick=false;
    bool AllQuad=false;
    bool Loaded=tri_mesh.LoadTriMesh(pathM,AllQuad);
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
    if (do_batch)
    {
        BatchProcess();
        exit(0);
    }
    //remeshed_mesh.UpdateDataStructures();
    Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
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
        tri_mesh.GLDrawSharpEdges();
        glWrap.Draw(vcg::GLW::DMFlatWire,vcg::GLW::CMNone,vcg::GLW::TMNone);
    }

    if (drawfield)
    {
        vcg::GLField<MyTriMesh>::GLDrawFaceField(tri_mesh,false,false,0.007);
        vcg::GLField<MyTriMesh>::GLDrawSingularity(tri_mesh);
    }

    if(hasToPick)
    {
        hasToPick=false;
        typename MyTriMesh::CoordType pp;
        if(vcg::Pick<typename MyTriMesh::CoordType>(xMouse,yMouse,pp))
        {
            typename MyTriMesh::CoordType closPt,bary;
            typename MyTriMesh::ScalarType minD;
            typename MyTriMesh::FaceType *f=vcg::tri::GetClosestFaceBase(tri_mesh,Gr,pp,tri_mesh.bbox.Diag(),minD,closPt);
            vcg::InterpolationParameters(*f,closPt,bary);
            size_t EdgeI=1;
            if ((bary.Y()<bary.X())&&(bary.Y()<bary.Z()))EdgeI=2;
            if ((bary.Z()<bary.X())&&(bary.Z()<bary.Y()))EdgeI=0;

            MyTriMesh::FaceType *fOpp=f->FFp(EdgeI);
            int eOpp=f->FFi(EdgeI);

            if (f->IsFaceEdgeS(EdgeI))
            {
                f->ClearFaceEdgeS(EdgeI);
                if (fOpp!=f)
                    fOpp->ClearFaceEdgeS(eOpp);
            }else
            {
                f->SetFaceEdgeS(EdgeI);
                if (fOpp!=f)
                    fOpp->SetFaceEdgeS(eOpp);
            }
        }
    }

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
