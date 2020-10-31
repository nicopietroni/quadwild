#include "generate_patch.h"

#include <patchgen/generate_topology.h>
#include <determine_geometry.h>

using namespace Eigen;

void patterns::generatePatch(const Eigen::VectorXi& l, patchgen::PatchParam& param, Patch& patch) {
    patchgen::generate_topology(l, param, patch);
    patterns::determine_geometry(patch, l);
}

void patterns::generatePatch(const patchgen::PatchParam& param, Patch& patch) {
    patchgen::generate_topology(param, patch);
    patterns::determine_geometry(patch, param.l);
}
