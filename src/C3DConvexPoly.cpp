#include "C3DConvexPoly.h"

#include <cmath>
#include <unordered_map>

namespace c3d {
namespace {

constexpr double kEpsilon = 1e-9;

Vec3 Normalize(const Vec3 &v) {
    double length = std::sqrt(Dot(v, v));
    if (length <= kEpsilon) {
        return Vec3{0.0, 0.0, 0.0};
    }
    return Vec3{v.x / length, v.y / length, v.z / length};
}

Face BuildFace(int a, int b, int c, const std::vector<Vec3> &vertices, const Vec3 &interior_point) {
    Face face;
    face.a = a;
    face.b = b;
    face.c = c;

    const Vec3 &va = vertices[a];
    const Vec3 &vb = vertices[b];
    const Vec3 &vc = vertices[c];

    Vec3 ab = vb - va;
    Vec3 ac = vc - va;
    face.normal = Normalize(Cross(ab, ac));
    face.offset = -Dot(face.normal, va);

    double interior_distance = Dot(face.normal, interior_point) + face.offset;
    if (interior_distance > 0.0) {
        std::swap(face.b, face.c);
        const Vec3 &vb2 = vertices[face.b];
        const Vec3 &vc2 = vertices[face.c];
        Vec3 ab2 = vb2 - va;
        Vec3 ac2 = vc2 - va;
        face.normal = Normalize(Cross(ab2, ac2));
        face.offset = -Dot(face.normal, va);
    }

    return face;
}

} // namespace

Face C3DConvexPoly::MakeFace(int a, int b, int c, const std::vector<Vec3> &vertices, const Vec3 &interior_point) {
    return BuildFace(a, b, c, vertices, interior_point);
}

C3DConvexPoly::C3DConvexPoly(double half_extent) {
    double h = half_extent;
    vertices_.reserve(16);
    vertices_.push_back(Vec3{-h, -h, -h});
    vertices_.push_back(Vec3{ h, -h, -h});
    vertices_.push_back(Vec3{ h,  h, -h});
    vertices_.push_back(Vec3{-h,  h, -h});
    vertices_.push_back(Vec3{-h, -h,  h});
    vertices_.push_back(Vec3{ h, -h,  h});
    vertices_.push_back(Vec3{ h,  h,  h});
    vertices_.push_back(Vec3{-h,  h,  h});

    interior_point_ = Vec3{0.0, 0.0, 0.0};

    faces_.reserve(12);
    auto add_face = [&](int a, int b, int c) {
        faces_.push_back(MakeFace(a, b, c, vertices_, interior_point_));
    };

    add_face(0, 1, 2);
    add_face(0, 2, 3);

    add_face(4, 6, 5);
    add_face(4, 7, 6);

    add_face(0, 4, 5);
    add_face(0, 5, 1);

    add_face(1, 5, 6);
    add_face(1, 6, 2);

    add_face(2, 6, 7);
    add_face(2, 7, 3);

    add_face(3, 7, 4);
    add_face(3, 4, 0);
}

bool C3DConvexPoly::Contains(const Vec3 &point, double epsilon) const {
    for (const Face &face : faces_) {
        double distance = Dot(face.normal, point) + face.offset;
        if (distance > epsilon) {
            return false;
        }
    }
    return true;
}

bool C3DConvexPoly::AddPoint(const Vec3 &point, double epsilon) {
    if (Contains(point, epsilon)) {
        return false;
    }

    const int new_index = static_cast<int>(vertices_.size());
    vertices_.push_back(point);
    double denom = static_cast<double>(new_index + 1);
    interior_point_ = (interior_point_ * static_cast<double>(new_index) + point) / denom;

    std::vector<int> visible_faces;
    visible_faces.reserve(faces_.size());
    for (int i = 0; i < static_cast<int>(faces_.size()); ++i) {
        const Face &face = faces_[i];
        double distance = Dot(face.normal, point) + face.offset;
        if (distance > epsilon) {
            visible_faces.push_back(i);
        }
    }

    if (visible_faces.empty()) {
        return false;
    }

    std::unordered_map<EdgeKey, int, EdgeKeyHasher> boundary_edges;
    auto add_edge = [&](int from, int to) {
        EdgeKey key{from, to};
        EdgeKey opposite{to, from};
        auto it = boundary_edges.find(opposite);
        if (it != boundary_edges.end()) {
            boundary_edges.erase(it);
        } else {
            boundary_edges.emplace(key, 1);
        }
    };

    std::vector<char> face_removed(faces_.size(), 0);
    for (int idx : visible_faces) {
        const Face &face = faces_[idx];
        add_edge(face.a, face.b);
        add_edge(face.b, face.c);
        add_edge(face.c, face.a);
        face_removed[idx] = 1;
    }

    std::vector<Face> new_faces;
    new_faces.reserve(boundary_edges.size());
    for (const auto &entry : boundary_edges) {
        int a = entry.first.from;
        int b = entry.first.to;
        Face face = MakeFace(a, b, new_index, vertices_, interior_point_);
        new_faces.push_back(face);
    }

    std::vector<Face> updated_faces;
    updated_faces.reserve(faces_.size() - visible_faces.size() + new_faces.size());
    for (int i = 0; i < static_cast<int>(faces_.size()); ++i) {
        if (!face_removed[i]) {
            updated_faces.push_back(faces_[i]);
        }
    }

    updated_faces.insert(updated_faces.end(), new_faces.begin(), new_faces.end());
    faces_.swap(updated_faces);

    return true;
}

} // namespace c3d

