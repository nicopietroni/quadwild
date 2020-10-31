#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <OpenMesh/Core/Mesh/Handles.hh>

// Modified and simplified version of <OpenMesh/Core/IO/reader/OBJReader.cc>

namespace kt84 {

template <class Mesh>
inline bool read_obj_soup(Mesh& mesh, const std::string& filename) {
    using namespace OpenMesh;

    std::ifstream fin( filename.c_str(), std::ios_base::in );
    if (!fin.is_open()) {
        std::cerr << "  Warning! Could not open file properly!\n";
        return false;
    }
    mesh.clear();
    
    std::string line;
    std::string keyWrd;

    std::vector<typename Mesh::Point> vertices;
    float                     x, y, z;

    auto trimString = [](std::string& _string) {
        // Trim Both leading and trailing spaces
        size_t start = _string.find_first_not_of(" \t\r\n");
        size_t end   = _string.find_last_not_of(" \t\r\n");
        if(( std::string::npos == start ) || ( std::string::npos == end))
            _string = "";
        else
            _string = _string.substr( start, end-start+1 );
    };
    
    while( fin && !fin.eof() )
    {
        std::getline(fin,line);
        if ( fin.bad() ){
            std::cerr << "  Warning! Could not read file properly!\n";
            return false;
        }

        // Trim Both leading and trailing spaces
        trimString(line);

        // comment
        if ( line.size() == 0 || line[0] == '#' || isspace(line[0]) ) {
            continue;
        }

        std::stringstream stream(line);

        stream >> keyWrd;

        // material file
        if (keyWrd == "mtllib")
        {
            // TODO
        }

        // usemtl
        else if (keyWrd == "usemtl")
        {
            // TODO
        }

        // vertex
        else if (keyWrd == "v")
        {
            stream >> x; stream >> y; stream >> z;

            if ( !stream.fail() )
            {
                vertices.push_back(typename Mesh::Point(x,y,z));
            }
        }

        // texture coord
        else if (keyWrd == "vt")
        {
            // TODO
        }

        // color per vertex
        else if (keyWrd == "vc")
        {
            // TODO
        }

        // normal
        else if (keyWrd == "vn")
        {
            // TODO
        }


        // face
        else if (keyWrd == "f")
        {
            int component(0), nV(0);
            int value;

            // read full line after detecting a face
            std::string faceLine;
            std::getline(stream,faceLine);
            std::stringstream lineData( faceLine );

            FaceHandle fh;
            std::vector<VertexHandle> faceVertices;

            // work on the line until nothing left to read
            while ( !lineData.eof() )
            {
                // read one block from the line ( vertex/texCoord/normal )
                std::string vertex;
                lineData >> vertex;

                do{

                    //get the component (vertex/texCoord/normal)
                    size_t found=vertex.find("/");

                    // parts are seperated by '/' So if no '/' found its the last component
                    if( found != std::string::npos ){

                        // read the index value
                        std::stringstream tmp( vertex.substr(0,found) );

                        // If we get an empty string this property is undefined in the file
                        if ( vertex.substr(0,found).empty() ) {
                            // Switch to next field
                            vertex = vertex.substr(found+1);

                            // Now we are at the next component
                            ++component;

                            // Skip further processing of this component
                            continue;
                        }

                        // Read current value
                        tmp >> value;

                        // remove the read part from the string
                        vertex = vertex.substr(found+1);

                    } else {

                        // last component of the vertex, read it.
                        std::stringstream tmp( vertex );
                        tmp >> value;

                        // Clear vertex after finished reading the line
                        vertex="";

                        // Nothing to read here ( garbage at end of line )
                        if ( tmp.fail() ) {
                            continue;
                        }
                    }

                    // store the component ( each component is referenced by the index here! )
                    switch (component)
                    {
                    case 0: // vertex
                        // Obj counts from 1 and not zero .. array counts from zero therefore -1
                        faceVertices.push_back(mesh.add_vertex(vertices[value-1]));
                        break;

                    case 1: // texture coord
                        // TODO
                        break;

                    case 2: // normal
                        // TODO
                        break;
                    }

                    // Prepare for reading next component
                    ++component;

                    // Read until line does not contain any other info
                } while ( !vertex.empty() );

                component = 0;
                nV++;

            }

            mesh.add_face(faceVertices);
        }
    }

    return true;
}

}
