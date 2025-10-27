#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace c3d {

// Simple value type representing a 3D vector.  The math helpers live next to the
// convex hull implementation because the hull only needs a subset of common
// vector operations.
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

// Triangular face of the hull (the hull is always stored as a triangle mesh so
// that it stays compatible with standard physics/graphics pipelines).
struct Face {
    int a;
    int b;
    int c;
    Vec3 normal;
    double offset;
};

// Incremental 3D convex hull maintainer.  The hull starts as a cube and each
// call to AddPoint() grows the hull only if the new point lies outside of the
// current volume.  The implementation follows the "incremental horizon"
// approach: visible faces are culled, the horizon edges are discovered, and new
// triangles are stitched between the horizon and the new vertex.
class C3DConvexPoly {
  public:
    explicit C3DConvexPoly(double half_extent = 1.0);

    // Inserts the point into the hull if it is outside the current volume.
    // Returns true only when a new vertex is accepted and the hull expands.
    bool AddPoint(const Vec3 &point, double epsilon = 1e-9);

    // Tests whether a point is inside (or on) the current hull within the
    // supplied tolerance.
    bool Contains(const Vec3 &point, double epsilon = 1e-9) const;

    // Rebuilds face normals and removes duplicate/degenerate faces that can be
    // created by floating-point noise or by importing a pre-existing model.
    // Calling this after a batch of AddPoint() calls keeps the hull minimal.
    void Optimize(double epsilon = 1e-9);

    const std::vector<Vec3> &Vertices() const { return vertices_; }
    const std::vector<Face> &Faces() const { return faces_; }

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
    static bool VerticesCollinear(const Vec3 &a, const Vec3 &b, const Vec3 &c, double epsilon);

    bool IsDuplicateVertex(const Vec3 &point, double epsilon) const;
    void PruneFaces(double epsilon);

    struct FaceKey {
        std::array<int, 3> indices;

        bool operator==(const FaceKey &other) const { return indices == other.indices; }
    };

    struct FaceKeyHasher {
        std::size_t operator()(const FaceKey &key) const noexcept {
            return static_cast<std::size_t>(key.indices[0]) * 73856093u ^
                   static_cast<std::size_t>(key.indices[1]) * 19349663u ^
                   static_cast<std::size_t>(key.indices[2]) * 83492791u;
        }
    };

    std::vector<Vec3> vertices_;
    std::vector<Face> faces_;
    Vec3 interior_point_;
};

} // namespace c3d

