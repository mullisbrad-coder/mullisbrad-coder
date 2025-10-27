#include "C3DConvexPoly.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>

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

double SquaredLength(const Vec3 &v) { return Dot(v, v); }

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

bool C3DConvexPoly::VerticesCollinear(const Vec3 &a, const Vec3 &b, const Vec3 &c, double epsilon) {
    Vec3 ab = b - a;
    Vec3 ac = c - a;
    Vec3 cross = Cross(ab, ac);
    return SquaredLength(cross) <= epsilon * epsilon;
}

bool C3DConvexPoly::IsDuplicateVertex(const Vec3 &point, double epsilon) const {
    double epsilon_sq = epsilon * epsilon;
    for (const Vec3 &vertex : vertices_) {
        Vec3 diff = vertex - point;
        if (SquaredLength(diff) <= epsilon_sq) {
            return true;
        }
    }
    return false;
}

void C3DConvexPoly::PruneFaces(double epsilon) {
    if (faces_.empty()) {
        return;
    }

    // Collapse duplicate index triplets and discard triangles that lost their
    // area because the input points became collinear.  This keeps the hull as
    // small as the requested tolerance allows.
    std::unordered_set<FaceKey, FaceKeyHasher> seen;
    seen.reserve(faces_.size());

    std::vector<Face> pruned;
    pruned.reserve(faces_.size());

    for (const Face &face : faces_) {
        if (face.a == face.b || face.b == face.c || face.c == face.a) {
            continue;
        }

        const Vec3 &va = vertices_[face.a];
        const Vec3 &vb = vertices_[face.b];
        const Vec3 &vc = vertices_[face.c];

        if (VerticesCollinear(va, vb, vc, epsilon)) {
            continue;
        }

        FaceKey key{std::array<int, 3>{face.a, face.b, face.c}};
        std::sort(key.indices.begin(), key.indices.end());
        if (!seen.emplace(key).second) {
            continue;
        }

        pruned.push_back(MakeFace(face.a, face.b, face.c, vertices_, interior_point_));
    }

    faces_.swap(pruned);
}

void C3DConvexPoly::CompactVertices() {
    if (vertices_.empty()) {
        interior_point_ = Vec3{0.0, 0.0, 0.0};
        return;
    }

    std::vector<char> used(vertices_.size(), 0);
    for (const Face &face : faces_) {
        used[face.a] = 1;
        used[face.b] = 1;
        used[face.c] = 1;
    }

    std::size_t active_count = 0;
    for (char flag : used) {
        active_count += static_cast<std::size_t>(flag);
    }

    if (active_count == 0) {
        vertices_.clear();
        interior_point_ = Vec3{0.0, 0.0, 0.0};
        return;
    }

    if (active_count == vertices_.size()) {
        Vec3 accum{0.0, 0.0, 0.0};
        for (const Vec3 &v : vertices_) {
            accum += v;
        }
        interior_point_ = accum / static_cast<double>(vertices_.size());
        return;
    }

    std::vector<Vec3> compact;
    compact.reserve(active_count);
    std::vector<int> remap(vertices_.size(), -1);

    for (std::size_t i = 0; i < vertices_.size(); ++i) {
        if (used[i]) {
            remap[i] = static_cast<int>(compact.size());
            compact.push_back(vertices_[i]);
        }
    }

    for (Face &face : faces_) {
        face.a = remap[face.a];
        face.b = remap[face.b];
        face.c = remap[face.c];
    }

    vertices_.swap(compact);

    Vec3 accum{0.0, 0.0, 0.0};
    for (const Vec3 &v : vertices_) {
        accum += v;
    }
    if (!vertices_.empty()) {
        interior_point_ = accum / static_cast<double>(vertices_.size());
    } else {
        interior_point_ = Vec3{0.0, 0.0, 0.0};
    }
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
    // Guard against inserting an equivalent vertex more than once—the horizon
    // algorithm assumes each vertex index maps to a unique position.
    if (IsDuplicateVertex(point, epsilon)) {
        return false;
    }

    if (Contains(point, epsilon)) {
        return false;
    }

    // Incremental horizon algorithm: remove all faces visible from the new
    // point, discover the border that separates the visible region from the
    // rest of the hull, and stitch new triangles between that border and the
    // vertex.
    const int new_index = static_cast<int>(vertices_.size());
    Vec3 previous_interior = interior_point_;
    vertices_.push_back(point);
    const double new_count = static_cast<double>(vertices_.size());
    interior_point_ = (previous_interior * (new_count - 1.0) + point) / new_count;

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
        // Numerical tolerances can disagree between Contains() and the face
        // visibility tests. When that happens we conservatively roll the
        // insertion back to keep the hull consistent and minimal.
        vertices_.pop_back();
        interior_point_ = previous_interior;
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
        // Stitch the horizon edge to the new vertex to form a counter-clockwise
        // triangle. MakeFace() orients the triangle so the normal still points
        // outwards from the hull.
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
    // Remove any degenerate or duplicate faces created by numerical error so
    // the hull stays as compact as possible.
    PruneFaces(epsilon);
    CompactVertices();

    return true;
}

void C3DConvexPoly::Optimize(double epsilon) {
    PruneFaces(epsilon);
    CompactVertices();
}

} // namespace c3d

