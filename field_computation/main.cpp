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

#include <QApplication>
#include <QDesktopWidget>
#include <GL/glew.h>
#include "glwidget.h"
#include <wrap/qt/anttweakbarMapper.h>
#include <QWindow>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QTextStream>
//#include <triangle_mesh_type.h>
//#include <wrap/io_trimesh/import_field.h>
//#include <wrap/io_trimesh/import.h>

#include<vcg/complex/algorithms/hole.h>


extern bool do_batch;
extern bool has_features;
extern bool has_features_fl;
extern std::string pathM,pathS;

extern bool do_remesh;
extern int remesher_iterations;
extern double remesher_aspect_ratio;
extern double alpha;
extern int remesher_termination_delta;

extern double sharp_feature_thr;
extern int feature_erode_dilate;

void correctOrAbort(const bool ok, const std::string & line, const int linenumber)
{
    if(!ok)
    {
        std::cerr << "[fieldComputation] Malformerd config file |-> " << line << " (" << linenumber << ")" << std::endl;
        std::cerr << "[fieldComputation] ...aborting..." << std::endl;
        abort();
    }
}

//bool loadConfigFile(const std::string & filename)
//{
//	QFile configFile(filename.c_str());

//	if (!configFile.open(QIODevice::ReadOnly))
//	{
//		//handleme
//	}

//	QTextStream configStream(&configFile);

//	int i = 0;
//    QString line;

//    while (configStream.readLineInto(&line))//(configStream.readLineInto(&line))
//	{
//		QStringList keyValuePair = line.split(' ');

//		if (keyValuePair.size() != 2)
//		{
//            correctOrAbort(false, line.toStdString(), i);
//		}

//		bool ok = false;
//		if (keyValuePair[0] == "remesh_iterations")
//		{
//			remesher_iterations = keyValuePair[1].toInt(&ok);
//            correctOrAbort(ok, line.toStdString(), i);
//			continue;
//		}
//		if (keyValuePair[0] == "remesh_target_aspect_ratio")
//		{
//			remesher_aspect_ratio = keyValuePair[1].toDouble(&ok);
//            correctOrAbort(ok, line.toStdString(), i);
//			continue;
//		}
//		if (keyValuePair[0] == "remesh_termination_delta")
//		{
//			remesher_termination_delta = keyValuePair[1].toInt(&ok);
//            correctOrAbort(ok, line.toStdString(), i);
//			continue;
//		}
//		if (keyValuePair[0] == "sharp_feature_thr")
//		{
//			sharp_feature_thr = keyValuePair[1].toDouble(&ok);
//            correctOrAbort(ok, line.toStdString(), i);
//			continue;
//		}
//		if (keyValuePair[0] == "sharp_feature_erode_dilate")
//		{
//			feature_erode_dilate = keyValuePair[1].toInt(&ok);
//            correctOrAbort(ok, line.toStdString(), i);
//			continue;
//		}

//		++i;
//	}

//	std::cout << "[fieldComputation] Successful config import" << std::endl;

//}

bool loadConfigFile(const std::string & filename)
{

    FILE *f=fopen(filename.c_str(),"rt");

    if (f==NULL)return false;

    std::cout<<"READ CONFIG FILE"<<std::endl;

    int IntVar;
    fscanf(f,"do_remesh %d\n",&IntVar);
    if (IntVar==0)
        do_remesh=false;
    std::cout<<"do_remesh "<<do_remesh<<std::endl;

    fscanf(f,"remesh_iterations %d\n",&remesher_iterations);
    std::cout<<"remesh_iterations "<<remesher_iterations<<std::endl;

    float remesher_aspect_ratiof;
    fscanf(f,"remesh_target_aspect_ratio %f\n",&remesher_aspect_ratiof);
    remesher_aspect_ratio=remesher_aspect_ratiof;
    std::cout<<"remesher_aspect_ratio "<<remesher_aspect_ratio<<std::endl;

    fscanf(f,"remesh_termination_delta %d\n",&remesher_termination_delta);
    std::cout<<"remesh_termination_delta "<<remesher_termination_delta<<std::endl;

    float sharp_feature_thrf;
    fscanf(f,"sharp_feature_thr %f\n",&sharp_feature_thrf);
    sharp_feature_thr=sharp_feature_thrf;
    std::cout<<"sharp_feature_thr "<<sharp_feature_thr<<std::endl;

    fscanf(f,"sharp_feature_erode_dilate %d\n",&feature_erode_dilate);
    std::cout<<"sharp_feature_erode_dilate "<<feature_erode_dilate<<std::endl;

    float alphaf;
    fscanf(f,"alpha %f\n",&alphaf);
    alpha=(double)alphaf;
    std::cout<<"alpha "<<alpha<<std::endl;

    fclose(f);

    std::cout << "[fieldComputation] Successful config import" << std::endl;

    return true;

}

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QWindow dummy;
    QString def_string = QString("GLOBAL fontscaling=%1").arg((int)dummy.devicePixelRatio());
    TwDefine(def_string.toStdString().c_str());
    printf("%s\n",qPrintable(def_string));
    fflush(stdout);

    // Set functions to handle string copy
    TwCopyCDStringToClientFunc(CopyCDStringToClient);
    TwCopyStdStringToClientFunc(CopyStdStringToClient);

    if( !TwInit(TW_OPENGL, NULL) )
    {
        fprintf(stderr, "AntTweakBar initialization failed: %s\n", TwGetLastError());
        return 1;
    }

    //load default config

    loadConfigFile("basic_setup.txt");
    assert(argc>1);
    pathM=std::string(argv[1]);
    for (size_t i=2;i<argc;i++)
    {
        if(std::string(argv[i])==std::string("batch"))
        {
            do_batch=true;
            continue;
        }

        std::string pathTest=std::string(argv[i]);
        int position=pathTest.find(".sharp");
        if (position!=-1)
        {
            pathS=pathTest;
            has_features=true;
            continue;
        }
        position=pathTest.find(".fl");
        if (position!=-1)
        {
            pathS=pathTest;
            has_features_fl=true;
            continue;
        }
        position=pathTest.find(".txt");
        if (position!=-1)
        {
           loadConfigFile(pathTest.c_str());
           continue;
        }
    }
    //    loadConfigFile("basic_setup.txt");

    //    std::cout << pathM << std::endl;

    //    if ((argc>2)&&(std::string(argv[2])==std::string("batch")))
    //        do_batch=true;
    //    else
    //    {
    //        if (argc>2)
    //        {
    //            std::string pathTest=std::string(argv[2]);
    //            int position=pathTest.find(".sharp");

    //            if (position!=-1)
    //            {
    //                pathS=pathTest;
    //                has_features=true;
    //            }
    //            position=pathTest.find(".fl");
    //            if (position!=-1)
    //            {
    //                pathS=pathTest;
    //                has_features_fl=true;
    //            }
    //        }
    //    }

    //    if ((has_features || has_features_fl)&&(argc>3))
    //    {
    //        if (std::string(argv[3])==std::string("batch"))
    //            do_batch=true;
    //    }
    GLWidget window;

    window.show();
    return app.exec();
}
