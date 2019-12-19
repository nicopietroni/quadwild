#ifndef MYUTILS_H
#define MYUTILS_H

#include "meshtypes.h"
#include <Eigen/Geometry>
#include <vcg/simplex/face/pos.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/refine.h>
#include <vcg/complex/algorithms/refine_loop.h>
#include <vcg/complex/algorithms/edge_collapse.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/inside.h>
#include <vcg/space/intersection2.h>
#include <vcg/space/triangle3.h>
#include <wrap/io_trimesh/export_ply.h>
#include <typeinfo>
#include <fstream>
namespace utility{

    //Verify is a pair of two vector are the same at least of a index permutation
    bool ispermuted(const vector<int>& list1,const vector<int>& list2, vector<pair<int,int>>& correspondence);

}

#endif // MYUTILS_H
