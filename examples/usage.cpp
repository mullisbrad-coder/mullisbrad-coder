#include "C3DConvexPoly.h"

#include <iostream>
#include <vector>

int main() {
    std::vector<c3d::Vec3> samples{
        {-1.5, -0.5, 0.25},
        {2.0, 0.0, 0.0},
        {0.0, 1.5, -0.5},
        {-0.75, -1.2, 1.0},
        {1.2, 0.8, 1.3},
        {-1.0, 1.1, -1.5},
        {0.5, -1.7, 0.8},
        {1.6, -0.6, -1.2}
    };

    c3d::C3DConvexPoly jarvis_hull = c3d::C3DConvexPoly::BuildJarvisMarch(samples);
    std::cout << "Jarvis hull vertices: " << jarvis_hull.Vertices().size() << "\n";
    std::cout << "Jarvis hull faces: " << jarvis_hull.Faces().size() << "\n";

    c3d::C3DConvexPoly graham_hull = c3d::C3DConvexPoly::BuildGrahamScan(samples);
    std::cout << "Graham hull vertices: " << graham_hull.Vertices().size() << "\n";
    std::cout << "Graham hull faces: " << graham_hull.Faces().size() << "\n";

    c3d::C3DConvexPoly dynamic_hull(1.0);
    for (const c3d::Vec3 &point : samples) {
        dynamic_hull.AddPoint(point);
    }

    std::cout << "Dynamic hull vertices: " << dynamic_hull.Vertices().size() << "\n";
    std::cout << "Dynamic hull faces: " << dynamic_hull.Faces().size() << "\n";
}

