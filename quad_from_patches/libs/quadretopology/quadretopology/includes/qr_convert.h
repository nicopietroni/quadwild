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

#ifndef QR_CONVERT_H
#define QR_CONVERT_H

#include <Eigen/Core>
#include <vector>

namespace QuadRetopology {
namespace internal {

template<class PolyMeshType>
void VCGToEigen(
        PolyMeshType& vcgMesh,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        std::vector<int>& vMap,
        std::vector<int>& fMap,
        bool selectedOnly = false,
        int numVerticesPerFace = 3,
        int dim = 3);

template<class PolyMeshType>
void eigenToVCG(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        PolyMeshType& vcgMesh,
        int numVertices = 3,
        int dim = 3);

}
}

#include "qr_convert.cpp"

#endif // QR_CONVERT_H
