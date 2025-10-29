#include "C3DConvexPoly.h"

#include <iostream>
#include <vector>

int main() {
    std::vector<c3d::Vec3> sample_points{
        {-2.0, -1.0, -0.5},
        {2.0, -1.0, -0.5},
        {2.0, 1.5, -0.5},
        {-2.0, 1.5, -0.5},
        {-2.0, -1.0, 0.5},
        {2.0, -1.0, 0.5},
        {2.0, 1.5, 0.5},
        {-2.0, 1.5, 0.5},
        {0.0, 0.0, 2.0},
        {0.5, -1.5, 1.0},
        {-1.0, 1.8, -1.5},
    };

    c3d::C3DConvexPoly dynamic_hull(1.0);
    for (const auto &point : sample_points) {
        dynamic_hull.AddPoint(point);
    }
    std::cout << "Dynamic incremental pre-optimisation vertices: " << dynamic_hull.Vertices().size()
              << " faces: " << dynamic_hull.Faces().size() << "\n";

    dynamic_hull.Optimize();

    auto jarvis_hull = c3d::BuildConvexPolyJarvis(sample_points);
    auto graham_hull = c3d::BuildConvexPolyGraham(sample_points);

    std::cout << "Dynamic incremental post-optimisation vertices: " << dynamic_hull.Vertices().size()
              << " faces: " << dynamic_hull.Faces().size() << "\n";
    std::cout << "Jarvis hull vertices: " << jarvis_hull.Vertices().size() << " faces: "
              << jarvis_hull.Faces().size() << "\n";
    std::cout << "Graham hull vertices: " << graham_hull.Vertices().size() << " faces: "
              << graham_hull.Faces().size() << "\n";
}

