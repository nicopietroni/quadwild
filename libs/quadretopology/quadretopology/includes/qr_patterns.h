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

#ifndef QR_PATTERNS_H
#define QR_PATTERNS_H

#include <vector>

#include <Eigen/Core>

namespace QuadRetopology {
namespace internal {

template<class PolyMesh>
void computePattern(
        const Eigen::VectorXi &l,
        Eigen::MatrixXd& patchV,
        Eigen::MatrixXi& patchF,
        PolyMesh& patchMesh,
        std::vector<size_t>& borders,
        std::vector<size_t>& corners,
        std::vector<std::vector<size_t>>& sides);

}
}

#include "qr_patterns.cpp"

#endif // QR_PATTERNS_H
