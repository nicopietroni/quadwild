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
#include <clocale>

#include "functions.h"

int main(int argc, char *argv[])
{
    Parameters parameters;

    FieldTriMesh trimesh;
    TraceMesh traceTrimesh;

    TriangleMesh trimeshToQuadrangulate;
    std::vector<std::vector<size_t>> trimeshPartitions;
    std::vector<std::vector<size_t>> trimeshCorners;
    std::vector<std::pair<size_t,size_t> > trimeshFeatures;
    std::vector<size_t> trimeshFeaturesC;

    PolyMesh quadmesh;
    std::vector<std::vector<size_t>> quadmeshPartitions;
    std::vector<std::vector<size_t>> quadmeshCorners;
    std::vector<int> ilpResult;

    if (argc==1)
    {
        std::cout<<"Please specify a mesh (OBJ or PLY) as argument"<<std::endl;
        exit(0);
    }

    std::string meshFilename=std::string(argv[1]);
    std::string sharpFilename;
    std::string fieldFilename;

    //Use "." as decimal separator
    std::setlocale(LC_NUMERIC, "en_US.UTF-8");

    std::cout<<"Reading input..."<<std::endl;
    loadConfigFile("basic_setup.txt", parameters);

    for (int i=2;i<argc;i++)
    {
        int position;

        std::string pathTest=std::string(argv[i]);

        position=pathTest.find(".sharp");
        if (position!=-1)
        {
            sharpFilename=pathTest;
            parameters.hasFeature=true;
            continue;
        }

        position=pathTest.find(".txt");
        if (position!=-1)
        {
           loadConfigFile(pathTest.c_str(), parameters);
           continue;
        }

        position=pathTest.find(".rosy");
        if (position!=-1)
        {
           fieldFilename=pathTest;
           parameters.hasField=true;
           continue;
        }
    }


    std::cout<<"Loading:"<<meshFilename.c_str()<<std::endl;

    bool allQuad;
    bool loaded=trimesh.LoadTriMesh(meshFilename,allQuad);
    trimesh.UpdateDataStructures();

    if (!loaded)
    {
        std::cout<<"Wrong mesh filename"<<std::endl;
        exit(0);
    }

    std::cout<<"Loaded "<<trimesh.fn<<" faces and "<<trimesh.vn<<" vertices"<<std::endl;

    std::cout<<std::endl<<"--------------------- 1 - Remesh and field ---------------------"<<std::endl;
    remeshAndField(trimesh, parameters, meshFilename, sharpFilename, fieldFilename);

    std::cout<<std::endl<<"--------------------- 2 - Tracing ---------------------"<<std::endl;
    trace(meshFilename, traceTrimesh);

    std::cout<<std::endl<<"--------------------- 3 - Quadrangulation ---------------------"<<std::endl;
    quadrangulate(meshFilename, trimeshToQuadrangulate, quadmesh, trimeshPartitions, trimeshCorners, trimeshFeatures, trimeshFeaturesC, quadmeshPartitions, quadmeshCorners, ilpResult, parameters);
}

