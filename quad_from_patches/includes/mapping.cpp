#include "mapping.h"

#include <igl/lscm.h>

#include <igl/AABB.h>
#include <igl/in_element.h>

#include <vcg/space/distance3.h>

namespace qfp {

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


void computeQuadrangulation(
        const Eigen::MatrixXd& chartV,
        const Eigen::MatrixXi& chartF,
        const Eigen::MatrixXd& patchV,
        const Eigen::MatrixXi& patchF,
        const std::vector<std::vector<size_t>>& chartSides,
        const std::vector<double>& chartSideLength,
        const std::vector<std::vector<size_t>>& patchSides,
        Eigen::MatrixXd& uvMapV,
        Eigen::MatrixXi& uvMapF,
        Eigen::MatrixXd& quadrangulationV,
        Eigen::MatrixXi& quadrangulationF)
{
    Eigen::VectorXi b;
    Eigen::MatrixXd bc;

    int chartBorderSize = 0;
    for (const std::vector<size_t>& side : chartSides)
        chartBorderSize += side.size()-1;

    b.resize(chartBorderSize);
    bc.resize(chartBorderSize, 2);

    int fixedId = 0;
    for (size_t sId = 0; sId < chartSides.size(); sId++) {
        //Get first and last corner of the side
        const size_t& firstPatchSideCornerId = patchSides[sId][0];
        const size_t& lastPatchSideCornerId = patchSides[sId][patchSides[sId].size() - 1];

        //Coordinate of the current corner
        const Eigen::VectorXd& cornerCoord = patchV.row(firstPatchSideCornerId);

        //Get vector of the side
        const Eigen::VectorXd vector = patchV.row(lastPatchSideCornerId) - patchV.row(firstPatchSideCornerId);
        double currentLength = 0;
        for (size_t i = 0; i < chartSides[sId].size() - 1; i++) {
            const std::vector<size_t>& chartSide = chartSides[sId];

            if (i > 0) {
                currentLength += (chartV.row(chartSide[i]) - chartV.row(chartSide[i-1])).norm();
            }

            size_t vId = chartSide[i];

            double lengthRatio = currentLength / chartSideLength[sId];
            assert(lengthRatio >= 0 && lengthRatio < 1);

            const Eigen::VectorXd uv = cornerCoord + (vector * lengthRatio);

            b(fixedId) = static_cast<int>(vId);

            //Flip x with y
            bc(fixedId, 0) = uv(1);
            bc(fixedId, 1) = uv(0);

            fixedId++;
        }
    }

    if (b.size() < chartV.rows()) {
        //Apply Least Square Conformal Maps
        igl::lscm(chartV, chartF, b, bc, uvMapV);
    }
    else {
        //Get the UV map with all fixed border
        uvMapV = bc;

#ifndef NDEBUG
        std::cout << "No fixed border! UVMap setted as bc." << std::endl;
#endif
    }

    uvMapF = chartF;

//    const Eigen::Vector2d& v1 = uvMapV.row(uvMapF(0,0));
//    const Eigen::Vector2d& v2 = uvMapV.row(uvMapF(0,1));
//    const Eigen::Vector2d& v3 = uvMapV.row(uvMapF(0,2));

//    Eigen::Vector3d e1(v2.x() - v1.x(), v2.y() - v1.y(), 0);
//    Eigen::Vector3d e2(v3.x() - v2.x(), v3.y() - v2.y(), 0);

//    Eigen::Vector3d normal = e1.cross(e2);

//    //Flip uv map faces
//    assert(normal.z() != 0);
//    if (normal.z() < 0) {
//        for (int i = 0; i < uvMapF.rows(); i++) {
//            assert(uvMapF.cols() == 3);
//            for (int j = 0; j < uvMapV.cols()/2; j++) {
//                std::swap(uvMapF(i,j), uvMapF(i, uvMapF.cols() - 1 - j));
//            }
//        }
//    }

    //AABB tree for point location
    igl::AABB<Eigen::MatrixXd, 2> tree;
    tree.init(uvMapV, uvMapF);

    quadrangulationV.resize(patchV.rows(), 3);
    for (int i = 0; i < patchV.rows(); i++) {
        const Eigen::VectorXd& Q = patchV.row(i);

        Eigen::MatrixXd Q2D(1,2);
        Q2D(0,0) = Q.x();
        Q2D(0,1) = Q.y();

        Eigen::VectorXi I;
        igl::in_element(uvMapV, uvMapF, Q2D, tree, I);

        int triIndex = I(0);

        if (triIndex < 0) {
            Eigen::VectorXi sqrD;
            Eigen::MatrixXd C;
            tree.squared_distance(uvMapV, uvMapF, Q2D, sqrD, I, C);

            triIndex = I(0);

            assert(triIndex > -1);
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

    //Flip faces
    for (int i = 0; i < quadrangulationF.rows(); i++) {
        assert(quadrangulationF.cols() == 4);
        for (int j = 0; j < quadrangulationF.cols()/2; j++) {
            std::swap(quadrangulationF(i,j), quadrangulationF(i, quadrangulationF.cols() - 1 - j));
        }
    }
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
        std::cout << "Degenerate triangle." << std::endl;
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

            bool wasCollinear = findParametricValueInSegment(t1,t2,closestPoint,paramT);
            assert(wasCollinear);

            baryc(0) = 1 - paramT.x();
            baryc(1) = paramT.x();
            baryc(2) = 0;
        }
        else if (dist2 <= dist1 && dist2 <= dist3) {
            closestPoint.x() = vcgClosestPoint2.X();
            closestPoint.y() = vcgClosestPoint2.Y();

            bool wasCollinear = findParametricValueInSegment(t2,t3,closestPoint,paramT);
            assert(wasCollinear);

            baryc(0) = 0;
            baryc(1) = 1 - paramT.x();
            baryc(2) = paramT.x();
        }
        else {
            closestPoint.x() = vcgClosestPoint3.X();
            closestPoint.y() = vcgClosestPoint3.Y();

            bool wasCollinear = findParametricValueInSegment(t3,t1,closestPoint,paramT);
            assert(wasCollinear);

            baryc(0) = paramT.x();
            baryc(1) = 0;
            baryc(2) = 1 - paramT.x();
        }

//        Eigen::VectorXd paramT(2);

//        if (findEquationInSegment(t1, t2, p, paramT)) {
//            baryc(0) = 1 - paramT.x();
//            baryc(1) = paramT.x();
//            baryc(2) = 0;
//#ifndef NDEBUG
//            std::cout << "Segment interpolation done.";
//#endif
//        }
//        else if (findEquationInSegment(t2, t3, p, paramT)) {
//            baryc(0) = 0;
//            baryc(1) = 1 - paramT.x();
//            baryc(2) = paramT.x();
//#ifndef NDEBUG
//            std::cout << "Segment interpolation done.";
//#endif
//        }
//        else if (findEquationInSegment(t3, t1, p, paramT)) {
//            baryc(0) = paramT.x();
//            baryc(1) = 0;
//            baryc(2) = 1 - paramT.x();
//#ifndef NDEBUG
//            std::cout << "Segment interpolation done.";
//#endif
//        }
//        else {
//            baryc(0) = baryc(1) = baryc(2) = 1.0/3.0;
//#ifndef NDEBUG
//            std::cout << "Segment interpolation NOT done." ;
//#endif
//        }
//#ifndef NDEBUG
//        std::cout << std::endl;
//#endif
    }
    else {
        baryc(2) = 1.0 - baryc(0) - baryc(1);
    }

    if (baryc(0) + baryc(1) + baryc(2) < 1.0 - eps || baryc(0) + baryc(1) + baryc(2) > 1.0 + eps) {
        baryc(0) = baryc(1) =  baryc(2) = 1.0/3.0;
        std::cout << "Barycenter sum not 1!" << std::endl;
    }

    if (std::isnan(baryc(0)) || std::isnan(baryc(1)) || std::isnan(baryc(2))) {
        baryc(0) = baryc(1) =  baryc(2) = 1.0/3.0;
        std::cout << "Barycenter nan!" << std::endl;
    }

    baryc(0) = std::max(std::min(baryc(0), 1.0), 0.0);
    baryc(1) = std::max(std::min(baryc(1), 1.0), 0.0);
    baryc(2) = std::max(std::min(baryc(2), 1.0), 0.0);

    assert(baryc(0) >= 0 && baryc(0) <= 1);
    assert(baryc(1) >= 0 && baryc(1) <= 1);
    assert(baryc(2) >= 0 && baryc(2) <= 1);


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

std::vector<std::vector<size_t>> getPatchSides(
        Eigen::MatrixXd& patchV,
        Eigen::MatrixXi& patchF,
        std::vector<size_t>& borders,
        std::vector<size_t>& corners,
        const Eigen::VectorXi& l)
{
    assert(corners.size() == l.size());

    std::vector<std::vector<size_t>> sides(corners.size());
    size_t startCornerId = 0;

    bool foundSolution;
    do {
        size_t bId = 0;
        while (borders[bId] != corners[startCornerId]) {
            bId = (bId + 1) % borders.size();
        }

        foundSolution = true;
        size_t cId = startCornerId;
        size_t sId = 0;
        do {
            assert(borders[bId] == corners[cId]);

            cId = (cId + 1) % corners.size();

            std::vector<size_t> side;

            while (borders[bId] != corners[cId]) {
                side.push_back(borders[bId]);
                bId = (bId + 1) % borders.size();
            }

            assert(borders[bId] == corners[cId]);
            assert(side[side.size() - 1] != corners[cId]);

            if (side.size() != l(sId)) {
                foundSolution = false;
#ifdef NDEBUG
                    break;
#endif
            }

            side.push_back(corners[cId]);

            sides[sId] = side;
            sId++;
        } while (startCornerId != cId);

       startCornerId++;
    } while (!foundSolution && startCornerId < corners.size());

    if (!foundSolution) {
#ifndef NDEBUG
      std::cout << "Found a no counter-clockwise path in patch. Reversed patch." << std::endl;
#endif

        for (int i = 0; i < patchV.rows(); i++) {
            for (int j = 0; j < patchV.cols(); j++) {
                patchV(i,j) = 0 - patchV(i,j);
            }
        }
        for (int i = 0; i < patchF.rows(); i++) {
            assert(patchF.cols() == 4);
            for (int j = 0; j < patchF.cols()/2; j++) {
                std::swap(patchF(i,j), patchF(i, patchF.cols() - 1 - j));
            }
        }

        std::reverse(corners.begin(), corners.end());
        std::reverse(borders.begin(), borders.end());

        startCornerId = 0;
        do {
            assert(startCornerId < corners.size());

            size_t bId = 0;
            while (borders[bId] != corners[startCornerId]) {
                bId = (bId + 1) % borders.size();
            }

            foundSolution = true;
            size_t cId = startCornerId;
            size_t sId = 0;
            do {
                assert(borders[bId] == corners[cId]);

                cId = (cId + 1) % corners.size();

                std::vector<size_t> side;

                while (borders[bId] != corners[cId]) {
                    side.push_back(borders[bId]);
                    bId = (bId + 1) % borders.size();
                }

                assert(borders[bId] == corners[cId]);
                assert(side[side.size() - 1] != corners[cId]);

                if (side.size() != l(sId)) {
                    foundSolution = false;
#ifdef NDEBUG
                    break;
#endif
                }

                side.push_back(corners[cId]);

                sides[sId] = side;
                sId++;
            } while (startCornerId != cId);

            startCornerId++;

        } while (!foundSolution);
    }

    assert(foundSolution);

    return sides;
}

}
