#pragma once
#include <patchgen/PatchParam.h>
#include "Patch.h"

namespace patterns {
    void generatePatch(const Eigen::VectorXi& l, patchgen::PatchParam& param, Patch& patch);
    void generatePatch(const patchgen::PatchParam& param, Patch& patch);
}
