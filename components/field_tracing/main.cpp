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

#include <QApplication>
#include <QDesktopWidget>
#include <wrap/qt/anttweakbarMapper.h>
#include "glwidget.h"
#include <QWindow>
#include <QFileInfo>

#include <clocale>

extern TraceMesh mesh;
extern std::string pathM;
extern std::string pathF;
extern std::string pathS;
extern std::string pathOF;
extern std::string pathProject;

extern bool has_features;
extern bool has_original_faces;
extern bool batch_process;

//extern float BatchSample;
//extern float BatchDrift;
//extern int BatchSplit;
//extern int BatchIncreaseValRem;
//extern int BatchIrregularRem;
//extern float BatchDistortionL;

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    //Use "." as decimal separator
    std::setlocale(LC_NUMERIC, "en_US.UTF-8");

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

    //PARAMETERS CHECK
    if(argc<2)
    {
        printf("error: pass mesh name as parameter \n");
        fflush(stdout);
        exit(0);
    }

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

    pathS=pathProject;
    pathS.append(".sharp");
    QString pathSQ=QString(pathS.c_str());
    QFileInfo f_infoS(pathSQ);
    if (!f_infoS.exists())
    {
        printf("no feature line \n");
        has_features=false;
        fflush(stdout);
        //exit(0);
    }
    else
    {
       has_features=true;
       std::cout<<"Sharp file correct"<<std::endl;
    }


    batch_process=false;
    if (argc>=3)
    {
        //then check if it mush batch process
        std::string pathComm;
        pathComm=std::string(argv[2]);
        if (pathComm==std::string("batch"))batch_process=true;
        if (batch_process)
        {
            std::cout<<"*** BATCH PROCESSING ***"<<std::endl;
            std::cout<<"* DATASET "<<pathM.c_str()<<" *"<<std::endl;
        }
    }


    pathOF=pathProject;
    pathOF.append("_origf.txt");

    QString pathOFQ=QString(pathOF.c_str());
    QFileInfo f_infoOF(pathOFQ);
    if (!f_infoOF.exists())
    {
        has_original_faces=false;
        printf("NO Original faces\n");
        fflush(stdout);
    }
    else
    {
        has_original_faces=true;
        std::cout<<"Original faces correct"<<std::endl;
    }

//    if ((batch_process)&&(argc>=4))
//        BatchSample=atof(argv[3]);

//    if ((batch_process)&&(argc>=5))
//        BatchDrift=atof(argv[4]);

//    if ((batch_process)&&(argc>=6))
//        BatchSplit=atoi(argv[5]);

//    if ((batch_process)&&(argc>=7))
//        BatchIncreaseValRem=atoi(argv[6]);

//    if ((batch_process)&&(argc>=8))
//        BatchIrregularRem=atoi(argv[7]);

//    if ((batch_process)&&(argc>=9))
//        BatchDistortionL=atof(argv[8]);

    GLWidget window;
    window.show();
    return app.exec();
}
