#include "qr_patterns.h"

#include "qr_convert.h"

#include <patterns/generate_patch.h>

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
        std::vector<std::vector<size_t>>& sides)
{
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
    patterns::Patch patch;
    patterns::generatePatch(l, param, patch);

    corners.resize(l.size());
    patchV.resize(patch.n_vertices(), 3);

    for (patterns::Patch::VertexIter v_it=patch.vertices_begin(); v_it!=patch.vertices_end(); ++v_it) {
        const int& corner_index = patch.data(*v_it).patchgen.corner_index;
        if (corner_index >= 0) {
            corners[corner_index] = v_it->idx();
        }

        patterns::Patch::Point p = patch.point(*v_it);

        patchV(v_it->idx(), 0) = p[0];
        patchV(v_it->idx(), 1) = p[1];
        patchV(v_it->idx(), 2) = p[2];
    }

    patchF.resize(patch.n_faces(), 4);
    for (patterns::Patch::FaceIter f_it=patch.faces_begin(); f_it!=patch.faces_end(); ++f_it) {
        int i = 0;
        for (patterns::Patch::FaceVertexIter fv_it = patch.fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
            patchF(f_it->idx(), i) = fv_it->idx();
            ++i;
        }
    }

    eigenToVCG(patchV, patchF, patchMesh, 4, 3);

    vcg::tri::UpdateTopology<PolyMesh>::FaceFace(patchMesh);

    std::vector<size_t> nextMap(patchV.rows(), patchV.rows() + 1);
    for (size_t i = 0; i < patchMesh.face.size(); i++) {
        for (int j = 0; j < patchMesh.face[i].VN(); j++) {
            if (vcg::face::IsBorder(patchMesh.face[i], j)) {
                size_t vStartId = vcg::tri::Index(patchMesh, patchMesh.face[i].V0(j));
                size_t vEndId = vcg::tri::Index(patchMesh, patchMesh.face[i].V1(j));
                nextMap[vStartId] = vEndId;
            }
        }
    }

    size_t currentId = corners[0];
    do {
        borders.push_back(currentId);
        currentId = nextMap[currentId];
    } while (currentId != corners[0]);


    size_t startCornerId = 0;

    size_t cId = 0;
    size_t bId = 0;
    assert(borders[bId] == corners[cId]);

    size_t sId = 0;
    sides.resize(corners.size());
    do {
        std::vector<size_t> side;

        cId = (cId + 1) % corners.size();

        while (borders[bId] != corners[cId]) {
            side.push_back(borders[bId]);
            bId = (bId + 1) % borders.size();
        }

        assert(borders[bId] == corners[cId]);
        assert(side[side.size() - 1] != corners[cId]);

        side.push_back(corners[cId]);

        sides[sId] = side;
        sId++;
    } while (cId != startCornerId);
}

}
}
