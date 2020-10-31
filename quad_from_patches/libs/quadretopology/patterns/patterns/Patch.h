#pragma once
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <patchgen/PatchVertexTraits.h>
#include <kt84/openmesh/base/LaplaceDirect.h>

namespace patterns {


struct PatchTraits : public OpenMesh::DefaultTraits {
        typedef OpenMesh::Vec3d Point;
        typedef OpenMesh::Vec3d Normal;

        VertexTraits
            , public patchgen::PatchVertexTraits
            , public kt84::LaplaceDirect_VertexTraits<3>
        {};

        HalfedgeTraits
            , public kt84::LaplaceDirect_HalfedgeTraits
        {};
};

typedef OpenMesh::PolyMesh_ArrayKernelT<PatchTraits> PatchBase;

struct Patch
    : public PatchBase
    , public kt84::LaplaceDirect<PatchBase, Patch, 3>
{
    /*
conventions for corner_index and boundary edge subdivision count l[i]

2-sided:
       l1
      ____
   __/    \__
  /          \
(0)          (1)
  \__      __/
     \____/
       l0

3-sided:
    (2)
    / \
l1 /   \  l2
  /     \
(0)-----(1)
    l0

4-sided:
      l2
 (3)------(2)
  |        |
l3|        |l1
  |        |
 (0)------(1)
      l0

5-sided
     __(3)__
l3 _/       \_ l2
  /           \
(4)           (2)
  \           /
l4 \         / l1
   (0)-----(1)
        l0

6-sided
         l3
    (4)-----(3)
    /         \
l4 /           \ l2
  /             \
(5)             (2)
  \             /
   \           / l1
l5  \         /
    (0)-----(1)
        l0
    */
};

}
