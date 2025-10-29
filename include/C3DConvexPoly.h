#pragma once

#include <cstddef>
#include <vector>

namespace c3d {

struct Vec3 {
    double x;
    double y;
    double z;

    Vec3() = default;
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 &operator+=(const Vec3 &rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    Vec3 &operator-=(const Vec3 &rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    Vec3 &operator*=(double s) {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    Vec3 operator+(const Vec3 &rhs) const { return Vec3{x + rhs.x, y + rhs.y, z + rhs.z}; }
    Vec3 operator-(const Vec3 &rhs) const { return Vec3{x - rhs.x, y - rhs.y, z - rhs.z}; }
    Vec3 operator*(double s) const { return Vec3{x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return Vec3{x / s, y / s, z / s}; }
};

inline Vec3 operator*(double s, const Vec3 &v) { return v * s; }

inline Vec3 Cross(const Vec3 &a, const Vec3 &b) {
    return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline double Dot(const Vec3 &a, const Vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

struct Face {
    int a;
    int b;
    int c;
    Vec3 normal;
    double offset;
};

class C3DConvexPoly {
  public:
    enum class HorizonTriangulation {
        kEdgeFan,
        kJarvisMarch,
        kGrahamScan,
    };

    explicit C3DConvexPoly(double half_extent = 1.0,
                           HorizonTriangulation method = HorizonTriangulation::kJarvisMarch);

    bool AddPoint(const Vec3 &point, double epsilon = 1e-9);

    bool Contains(const Vec3 &point, double epsilon = 1e-9) const;

    const std::vector<Vec3> &Vertices() const { return vertices_; }
    const std::vector<Face> &Faces() const { return faces_; }

    void SetTriangulationMethod(HorizonTriangulation method) { method_ = method; }
    HorizonTriangulation TriangulationMethod() const { return method_; }

  private:
    struct EdgeKey {
        int from;
        int to;

        bool operator==(const EdgeKey &other) const { return from == other.from && to == other.to; }
    };

    struct EdgeKeyHasher {
        std::size_t operator()(const EdgeKey &edge) const noexcept {
            return static_cast<std::size_t>(edge.from) * 1315423911u + static_cast<std::size_t>(edge.to);
        }
    };

    static Face MakeFace(int a, int b, int c, const std::vector<Vec3> &vertices, const Vec3 &interior_point);

    struct ProjectedVertex {
        int index;
        double x;
        double y;
    };

    std::vector<Face> BuildHorizonFaces(int new_index, const std::vector<EdgeKey> &edges) const;
    std::vector<int> JarvisMarchHull(const std::vector<ProjectedVertex> &points) const;
    std::vector<int> GrahamScanHull(const std::vector<ProjectedVertex> &points) const;
    std::vector<ProjectedVertex> ProjectComponent(const std::vector<int> &component_vertices,
                                                  const std::vector<int> &component_edges,
                                                  int new_index,
                                                  const std::vector<EdgeKey> &edges) const;
    static double Orientation2D(const ProjectedVertex &a, const ProjectedVertex &b, const ProjectedVertex &c);
    static double DistanceSq2D(const ProjectedVertex &a, const ProjectedVertex &b);

    std::vector<Vec3> vertices_;
    std::vector<Face> faces_;
    Vec3 interior_point_;
    HorizonTriangulation method_;
};

} // namespace c3d

