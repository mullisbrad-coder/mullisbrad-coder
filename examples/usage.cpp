#include "C3DConvexPoly.h"

#include <iostream>

int main() {
    c3d::C3DConvexPoly hull(1.0);

    c3d::Vec3 p1{2.0, 0.0, 0.0};
    c3d::Vec3 p2{0.0, 0.5, 0.0};

    bool added1 = hull.AddPoint(p1);
    bool added2 = hull.AddPoint(p2);

    std::cout << "Added1: " << added1 << "\n";
    std::cout << "Added2: " << added2 << "\n";
    std::cout << "Total vertices: " << hull.Vertices().size() << "\n";
    std::cout << "Total faces: " << hull.Faces().size() << "\n";
}

