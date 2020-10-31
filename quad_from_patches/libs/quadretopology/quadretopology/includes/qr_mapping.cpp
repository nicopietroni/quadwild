#include "qr_mapping.h"

#include <igl/lscm.h>

#include <igl/AABB.h>
#include <igl/in_element.h>

#include <vcg/space/distance3.h>

namespace QuadRetopology {
namespace internal {

Eigen::VectorXd pointToBarycentric(
        const Eigen::VectorXd& t1,
        const Eigen::VectorXd& t2,
        const Eigen::VectorXd& t3,
        const Eigen::VectorXd& p);

Eigen::VectorXd barycentricToPoint(
        const Eigen::VectorXd& t1,
        const Eigen::VectorXd& t2,
        const Eigen::VectorXd& t3,
        const Eigen::VectorXd& p);

bool findParametricValueInSegment(
        const Eigen::VectorXd& s1,
        const Eigen::VectorXd& s2,
        const Eigen::VectorXd& p,
        Eigen::VectorXd& t);

void computeQuadrangulation(
        const Eigen::MatrixXd& chartV,
        const Eigen::MatrixXi& chartF,
        const Eigen::MatrixXd& patchV,
        const Eigen::MatrixXi& patchF,
        const std::vector<std::vector<std::vector<size_t>>>& chartSideVertices,
        const std::vector<std::vector<double>>& chartSideLength,
        const std::vector<std::vector<size_t>>& chartSideSubdivision,
        const std::vector<std::vector<size_t>>& patchSides,
        Eigen::MatrixXd& uvMapV,
        Eigen::MatrixXi& uvMapF,
        Eigen::MatrixXd& quadrangulationV,
        Eigen::MatrixXi& quadrangulationF)
{
    Eigen::VectorXi b;
    Eigen::MatrixXd bc;

    int chartBorderSize = 0;
    for (const std::vector<std::vector<size_t>>& side : chartSideVertices) {
        for (const std::vector<size_t>& subside : side) {
            chartBorderSize += subside.size()-1;
        }
    }

    b.resize(chartBorderSize);
    bc.resize(chartBorderSize, 2);

    int fixedId = 0;
    for (size_t i = 0; i < chartSideVertices.size(); i++) {
        size_t patchStart = 0;

        for (size_t j = 0; j < chartSideVertices[i].size(); j++) {
            size_t patchEnd = patchStart + chartSideSubdivision[i][j];

            //Get first and last corner of the side
            const size_t& patchFirstId = patchSides[i][patchStart];
            const size_t& patchLastId = patchSides[i][patchEnd];

            //Coordinate of the current corner
            const Eigen::VectorXd& patchFirstCoord = patchV.row(patchFirstId);
            const Eigen::VectorXd& patchLastCoord = patchV.row(patchLastId);

            //Get vector of the side
            const Eigen::VectorXd vector = patchLastCoord - patchFirstCoord;

            double currentLength = 0;
            for (size_t k = 0; k < chartSideVertices[i][j].size() - 1; k++) {
                if (k > 0) {
                    currentLength += (chartV.row(chartSideVertices[i][j][k]) - chartV.row(chartSideVertices[i][j][k-1])).norm();
                }

                size_t vId = chartSideVertices[i][j][k];

                double lengthRatio = currentLength / chartSideLength[i][j];
                assert(lengthRatio >= 0 && lengthRatio < 1);

                const Eigen::VectorXd uv = patchFirstCoord + (vector * lengthRatio);

                b(fixedId) = static_cast<int>(vId);

                bc(fixedId, 0) = uv(0);
                bc(fixedId, 1) = uv(1);

                fixedId++;
            }

            patchStart += chartSideSubdivision[i][j];
        }
    }

    if (b.size() < chartV.rows()) {
        //Apply Least Square Conformal Maps
        igl::lscm(chartV, chartF, b, bc, uvMapV);

        //Flip x with y (problem with libigl)
        for (int i = 0; i < uvMapV.rows(); i++) {
            uvMapV.row(i).reverseInPlace();
        }
    }
    else {
        //Get the UV map with all fixed border
        uvMapV.resize(chartV.rows(), 2);
        for (int i = 0; i < bc.rows(); i++) {
            uvMapV.row(b(i)) = bc.row(i);
        }

#ifndef NDEBUG
        std::cout << "No fixed border! UVMap setted as bc." << std::endl;
#endif
    }

    uvMapF = chartF;

#ifndef NDEBUG
    for (int i = 0; i < uvMapF.rows(); ++i) {
        Eigen::VectorXd v1 = uvMapV.row(uvMapF(i,0));
        Eigen::VectorXd v2 = uvMapV.row(uvMapF(i,1));
        Eigen::VectorXd v3 = uvMapV.row(uvMapF(i,2));
        double A = std::abs(
            (v1.x()*v2.y() - v2.x()*v1.y()) +
            (v2.x()*v3.y() - v3.x()*v2.y()) +
            (v3.x()*v1.y() - v1.x()*v3.y())
        ) / 2.0;
        if (A <= 0) {
            std::cout << "Warning: degenerate triangle found in UV mapping!" << std::endl;
        }
    }
#endif

    //AABB tree for point location
    igl::AABB<Eigen::MatrixXd, 2> tree;
    tree.init(uvMapV, uvMapF);

    quadrangulationV.resize(patchV.rows(), 3);
    for (int i = 0; i < patchV.rows(); i++) {
        const Eigen::VectorXd& Q = patchV.row(i);

        Eigen::MatrixXd Q2D(1,2);
        Q2D(0,0) = Q.x();
        Q2D(0,1) = Q.y();

        int triIndex;
        if (uvMapF.rows() == 1) {
            triIndex = 0;
        }
        else {
            Eigen::VectorXi I;
            igl::in_element(uvMapV, uvMapF, Q2D, tree, I);

            triIndex = I(0);

            if (triIndex < 0) {
                Eigen::VectorXd sqrD;
                Eigen::MatrixXd C;
                tree.squared_distance(uvMapV, uvMapF, Q2D, sqrD, I, C);

                double currentSquareDistance = sqrD(0);
                triIndex = I(0);
                for (int k = 1; k < I.size(); k++) {
                    if (sqrD(k) < currentSquareDistance) {
                        triIndex = I(k);
                        currentSquareDistance = sqrD(k);
                    }
                }

                assert(triIndex > -1);
            }
        }

        const Eigen::VectorXi& tri = chartF.row(triIndex);

        Eigen::VectorXd baryc = pointToBarycentric(
                    uvMapV.row(tri(0)),
                    uvMapV.row(tri(1)),
                    uvMapV.row(tri(2)),
                    Q2D.row(0));

        Eigen::VectorXd mappedPoint = barycentricToPoint(
                    chartV.row(tri(0)),
                    chartV.row(tri(1)),
                    chartV.row(tri(2)),
                    baryc);

        quadrangulationV.row(i) = mappedPoint;
    }

    quadrangulationF = patchF;
}


Eigen::VectorXd pointToBarycentric(
        const Eigen::VectorXd& t1,
        const Eigen::VectorXd& t2,
        const Eigen::VectorXd& t3,
        const Eigen::VectorXd& p)
{
    const double eps = 0.0001;
    double det = (t2.y() - t3.y()) * (t1.x() - t3.x()) + (t3.x() - t2.x()) * (t1.y() - t3.y());

    Eigen::VectorXd baryc(3);

    baryc(0) = ((t2.y() - t3.y()) * (p.x() - t3.x()) + (t3.x() - t2.x()) * (p.y() - t3.y())) / det;
    baryc(1) = ((t3.y() - t1.y()) * (p.x() - t3.x()) + (t1.x() - t3.x()) * (p.y() - t3.y())) / det;

    if (baryc(0) > 1.0 + eps || baryc(1) > 1.0 + eps || baryc(0) < 0.0 - eps || baryc(1) < 0.0 - eps) {
#ifndef NDEBUG
        std::cout << "Valid barycenter coordinates not found." << std::endl;
#endif

        vcg::Point3d t1vcg(t1.x(), t1.y(), 0);
        vcg::Point3d t2vcg(t2.x(), t2.y(), 0);
        vcg::Point3d t3vcg(t3.x(), t3.y(), 0);
        vcg::Point3d pvcg(p.x(), p.y(), 0);
        vcg::Segment3d seg1(t1vcg, t2vcg);
        vcg::Segment3d seg2(t2vcg, t3vcg);
        vcg::Segment3d seg3(t3vcg, t1vcg);

        vcg::Point3d vcgClosestPoint1;
        vcg::Point3d vcgClosestPoint2;
        vcg::Point3d vcgClosestPoint3;
        double dist1;
        double dist2;
        double dist3;

        vcg::SegmentPointDistance(seg1, pvcg, vcgClosestPoint1, dist1);
        vcg::SegmentPointDistance(seg2, pvcg, vcgClosestPoint2, dist2);
        vcg::SegmentPointDistance(seg3, pvcg, vcgClosestPoint3, dist3);

        Eigen::VectorXd paramT(2);
        Eigen::VectorXd closestPoint(2);
        if (dist1 <= dist2 && dist1 <= dist3) {
            closestPoint.x() = vcgClosestPoint1.X();
            closestPoint.y() = vcgClosestPoint1.Y();

            findParametricValueInSegment(t1,t2,closestPoint,paramT);

            baryc(0) = 1 - paramT.x();
            baryc(1) = paramT.x();
            baryc(2) = 0;
        }
        else if (dist2 <= dist1 && dist2 <= dist3) {
            closestPoint.x() = vcgClosestPoint2.X();
            closestPoint.y() = vcgClosestPoint2.Y();

            findParametricValueInSegment(t2,t3,closestPoint,paramT);

            baryc(0) = 0;
            baryc(1) = 1 - paramT.x();
            baryc(2) = paramT.x();
        }
        else {
            closestPoint.x() = vcgClosestPoint3.X();
            closestPoint.y() = vcgClosestPoint3.Y();

            findParametricValueInSegment(t3,t1,closestPoint,paramT);

            baryc(0) = paramT.x();
            baryc(1) = 0;
            baryc(2) = 1 - paramT.x();
        }
    }
    else {
        baryc(2) = 1.0 - baryc(0) - baryc(1);
    }

    if (baryc(0) + baryc(1) + baryc(2) < 1.0 - eps || baryc(0) + baryc(1) + baryc(2) > 1.0 + eps) {
        baryc(0) = baryc(1) =  baryc(2) = 1.0/3.0;
#ifndef NDEBUG
        std::cout << "Barycenter sum not 1!" << std::endl;
#endif
    }

    if (std::isnan(baryc(0)) || std::isnan(baryc(1)) || std::isnan(baryc(2))) {
        baryc(0) = baryc(1) =  baryc(2) = 1.0/3.0;
#ifndef NDEBUG
        std::cout << "Barycenter nan!" << std::endl;
#endif
    }

    baryc(0) = std::max(std::min(baryc(0), 1.0), 0.0);
    baryc(1) = std::max(std::min(baryc(1), 1.0), 0.0);
    baryc(2) = std::max(std::min(baryc(2), 1.0), 0.0);

    return baryc;
}

Eigen::VectorXd barycentricToPoint(
        const Eigen::VectorXd& t1,
        const Eigen::VectorXd& t2,
        const Eigen::VectorXd& t3,
        const Eigen::VectorXd& baryc)
{
    Eigen::VectorXd coordinates =
            t1 * baryc(0) +
            t2 * baryc(1) +
            t3 * baryc(2);

    return coordinates;
}

bool findParametricValueInSegment(
        const Eigen::VectorXd& s1,
        const Eigen::VectorXd& s2,
        const Eigen::VectorXd& p,
        Eigen::VectorXd& t)
{
    const double eps = 0.0001;

    t.x() = (p.x() - s1.x())/(s2.x() - s1.x());
    t.y() = (p.y() - s1.y())/(s2.y() - s1.y());

    //Return true if it is collinear and inside the segment
    return  t.x() >= 0 - eps &&
            t.x() <= 1 + eps &&
            t.y() >= 0 - eps &&
            t.y() <= 1 + eps &&
            t.x() - eps <= t.y() + eps &&
            t.x() + eps >= t.y() - eps;

}

}
}
