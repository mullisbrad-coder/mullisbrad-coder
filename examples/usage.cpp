#include "C3DConvexPoly.h"

#include <iostream>
#include <string>
#include <vector>

void PrintHull(const std::string &label, const c3d::C3DConvexPoly &hull) {
    std::cout << label << "\n";
    std::cout << "  Vertices: " << hull.Vertices().size() << "\n";
    std::cout << "  Faces: " << hull.Faces().size() << "\n";
}

int main() {
    using c3d::C3DConvexPoly;
    using c3d::Vec3;

    std::vector<Vec3> samples = {
        Vec3{2.0, 0.0, 0.0},
        Vec3{0.0, 2.0, 0.0},
        Vec3{0.0, 0.0, 2.0},
        Vec3{-1.5, -1.0, 1.2},
        Vec3{1.0, -1.5, -1.1},
    };

    C3DConvexPoly jarvis_hull(1.0, C3DConvexPoly::HorizonTriangulation::kJarvisMarch);
    C3DConvexPoly graham_hull(1.0, C3DConvexPoly::HorizonTriangulation::kGrahamScan);

    for (const Vec3 &p : samples) {
        jarvis_hull.AddPoint(p);
        graham_hull.AddPoint(p);
    }

    PrintHull("Jarvis March hull", jarvis_hull);
    PrintHull("Graham Scan hull", graham_hull);

    return 0;
}
