#include "decl.h"
using namespace std;
using namespace Eigen;

Vector2d patchgen::get_boundary_geometry(int num_sides, double t) {
    const double pi = 3.1415926535897932384626433832795;

    if (num_sides == 1) {
        double theta = t * 2 * pi;
        return Vector2d(cos(theta), sin(theta));            // I don't know what shape is expected here.
    }

    if (num_sides == 2) {
        if (t < 1) {
            // left to right
            double theta = (2 * t - 1) * 3.14159 / 3;
            return Vector2d(sin(theta), -cos(theta) + 0.5);
        } else {
            t -= 1;
            // right to left
            double theta = (2 * t - 1) * 3.14159 / 3;
            return Vector2d(-sin(theta), cos(theta) - 0.5);
        }
    }

    double side_angle = 2 * pi / num_sides;
    vector<Vector2d> corners(num_sides);
    for (int i = 0; i < num_sides; ++i) {
        double theta = side_angle * i;
        // orient such that the first corner comes to the bottom left
        theta -= 0.5 * pi + 0.5 * side_angle;
        corners[i] << cos(theta), sin(theta);
    }

    if (t == num_sides) return corners[0];

    int i1 = static_cast<int>(t);
    int i2 = (i1 + 1) % num_sides;
    double s = t - i1;
    return (1 - s) * corners[i1] + s * corners[i2];
}
