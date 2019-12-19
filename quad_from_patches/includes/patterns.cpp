#include "patterns.h"

#include <patterns/ktmethod/patchgen/PatchParam.h>
#include <patterns/meshtypes.h>
#include <patterns/patchg.h>

#include <iostream>

#include "convert.h"

namespace qfp {

void computePattern(
        const Eigen::VectorXi &l,
        Eigen::MatrixXd& patchV,
        Eigen::MatrixXi& patchF,
        std::vector<size_t>& borders,
        std::vector<size_t>& corners)
{
    PatchG<PMesh,patchgen::PatchParam> patch;

    long int num_sides = l.size();
    if (num_sides < 2 || 6 < num_sides) {
        std::cout << "num_sides=" << num_sides << " is unsupported.\n";
        return;
    }
    if (l.sum() % 2 != 0) {
        std::cout << "The sum of number of edge subdivisions should be even.\n";
        return;
    }
    if (l.sum() < 4) {
        std::cout << "Input numbers are too small.\n";
        return;
    }

    patchgen::PatchParam param;
    patch.generate_topology(l, param);
    patch.determine_geometry(l);

    std::vector<int> vMap;
    std::vector<int> fMap;
    VCGToEigen(patch.mesh, patchV, patchF, vMap, fMap, false, 4);

    for (int& borderId : patch.borders) {
        borders.push_back(vMap[borderId]);
    }
    for (int& cornerId : patch.corners) {
        corners.push_back(vMap[cornerId]);
    }
}

}
