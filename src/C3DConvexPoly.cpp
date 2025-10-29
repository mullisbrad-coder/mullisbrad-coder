#include "C3DConvexPoly.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace c3d {
namespace {

constexpr double kEpsilon = 1e-9;
constexpr double kTwoPi = 6.283185307179586476925286766559;

struct FaceKey {
    int a;
    int b;
    int c;

    bool operator==(const FaceKey &other) const {
        return a == other.a && b == other.b && c == other.c;
    }
};

struct FaceKeyHasher {
    std::size_t operator()(const FaceKey &key) const noexcept {
        std::size_t seed = static_cast<std::size_t>(key.a);
        seed = seed * 1315423911u + static_cast<std::size_t>(key.b);
        seed = seed * 1315423911u + static_cast<std::size_t>(key.c);
        return seed;
    }
};

struct DirectedEdge {
    int from;
    int to;
    int previous;
};

struct DirectedEdgeKey {
    int from;
    int to;

    bool operator==(const DirectedEdgeKey &other) const {
        return from == other.from && to == other.to;
    }
};

struct DirectedEdgeHasher {
    std::size_t operator()(const DirectedEdgeKey &edge) const noexcept {
        return static_cast<std::size_t>(edge.from) * 1315423911u + static_cast<std::size_t>(edge.to);
    }
};

struct TetrahedronData {
    bool valid{false};
    std::array<int, 4> indices{0, 1, 2, 3};
    std::vector<Vec3> vertices;
    std::vector<Face> faces;
    Vec3 interior_point{0.0, 0.0, 0.0};
};

double LengthSquared(const Vec3 &v) {
    return Dot(v, v);
}

Vec3 Normalize(const Vec3 &v) {
    double length = std::sqrt(LengthSquared(v));
    if (length <= kEpsilon) {
        return Vec3{0.0, 0.0, 0.0};
    }
    return Vec3{v.x / length, v.y / length, v.z / length};
}

double SignedVolume(const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &d) {
    return Dot(Cross(b - a, c - a), d - a);
}

double DistanceToLineSquared(const Vec3 &a, const Vec3 &b, const Vec3 &p) {
    Vec3 ab = b - a;
    Vec3 ap = p - a;
    Vec3 cross = Cross(ab, ap);
    double denom = LengthSquared(ab);
    if (denom <= kEpsilon) {
        return LengthSquared(ap);
    }
    return LengthSquared(cross) / denom;
}

FaceKey CanonicalFaceKey(int a, int b, int c) {
    std::array<int, 3> values{a, b, c};
    std::sort(values.begin(), values.end());
    return FaceKey{values[0], values[1], values[2]};
}

bool IsFaceDegenerate(const Face &face) {
    return LengthSquared(face.normal) <= kEpsilon;
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

TetrahedronData BuildInitialTetrahedron(const std::vector<Vec3> &points, double epsilon) {
    TetrahedronData data;
    if (points.size() < 4) {
        return data;
    }

    const int count = static_cast<int>(points.size());
    int p0 = 0;
    for (int i = 1; i < count; ++i) {
        if (points[i].x < points[p0].x ||
            (std::abs(points[i].x - points[p0].x) <= epsilon && points[i].y < points[p0].y) ||
            (std::abs(points[i].x - points[p0].x) <= epsilon && std::abs(points[i].y - points[p0].y) <= epsilon &&
             points[i].z < points[p0].z)) {
            p0 = i;
        }
    }

    int p1 = -1;
    double best_dist = -1.0;
    for (int i = 0; i < count; ++i) {
        if (i == p0) {
            continue;
        }
        double dist = LengthSquared(points[i] - points[p0]);
        if (dist > best_dist) {
            best_dist = dist;
            p1 = i;
        }
    }

    if (p1 == -1 || best_dist <= epsilon) {
        return data;
    }

    int p2 = -1;
    double best_area = -1.0;
    for (int i = 0; i < count; ++i) {
        if (i == p0 || i == p1) {
            continue;
        }
        Vec3 normal = Cross(points[p1] - points[p0], points[i] - points[p0]);
        double area = std::sqrt(LengthSquared(normal));
        if (area > best_area) {
            best_area = area;
            p2 = i;
        }
    }

    if (p2 == -1 || best_area <= epsilon) {
        return data;
    }

    int p3 = -1;
    double best_volume = -1.0;
    for (int i = 0; i < count; ++i) {
        if (i == p0 || i == p1 || i == p2) {
            continue;
        }
        double volume = std::abs(SignedVolume(points[p0], points[p1], points[p2], points[i]));
        if (volume > best_volume) {
            best_volume = volume;
            p3 = i;
        }
    }

    if (p3 == -1 || best_volume <= epsilon) {
        return data;
    }

    data.valid = true;
    data.indices = {p0, p1, p2, p3};
    data.vertices = {points[p0], points[p1], points[p2], points[p3]};
    data.interior_point = (points[p0] + points[p1] + points[p2] + points[p3]) * 0.25;

    data.faces.reserve(4);
    Face f0 = BuildFace(0, 1, 2, data.vertices, data.interior_point);
    Face f1 = BuildFace(0, 3, 1, data.vertices, data.interior_point);
    Face f2 = BuildFace(1, 3, 2, data.vertices, data.interior_point);
    Face f3 = BuildFace(2, 3, 0, data.vertices, data.interior_point);

    if (IsFaceDegenerate(f0) || IsFaceDegenerate(f1) || IsFaceDegenerate(f2) || IsFaceDegenerate(f3)) {
        data.valid = false;
        data.faces.clear();
        data.vertices.clear();
        return data;
    }

    data.faces.push_back(f0);
    data.faces.push_back(f1);
    data.faces.push_back(f2);
    data.faces.push_back(f3);

    return data;
}

int FindGiftWrapCandidate(int from, int to, int previous, const std::vector<Vec3> &points, double epsilon) {
    const Vec3 &pa = points[from];
    const Vec3 &pb = points[to];
    Vec3 edge = pb - pa;

    int candidate = -1;
    double best_distance = -1.0;
    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        if (i == from || i == to || i == previous) {
            continue;
        }

        if (candidate == -1) {
            candidate = i;
            best_distance = DistanceToLineSquared(pa, pb, points[i]);
            continue;
        }

        double volume = SignedVolume(pa, pb, points[candidate], points[i]);
        if (volume < -epsilon) {
            candidate = i;
            best_distance = DistanceToLineSquared(pa, pb, points[i]);
        } else if (std::abs(volume) <= epsilon) {
            double current_distance = DistanceToLineSquared(pa, pb, points[i]);
            if (current_distance > best_distance + kEpsilon) {
                candidate = i;
                best_distance = current_distance;
            }
        }
    }

    return candidate;
}

} // namespace

Face C3DConvexPoly::MakeFace(int a, int b, int c, const std::vector<Vec3> &vertices, const Vec3 &interior_point) {
    return BuildFace(a, b, c, vertices, interior_point);
}

C3DConvexPoly::C3DConvexPoly(std::vector<Vec3> vertices, std::vector<Face> faces, Vec3 interior_point)
    : vertices_(std::move(vertices)), faces_(std::move(faces)), interior_point_(interior_point) {}

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

    std::unordered_set<FaceKey, FaceKeyHasher> face_keys;
    face_keys.reserve(faces_.size());
    for (const Face &face : faces_) {
        face_keys.insert(CanonicalFaceKey(face.a, face.b, face.c));
    }

    std::vector<char> face_removed(faces_.size(), 0);
    for (int idx : visible_faces) {
        const Face &face = faces_[idx];
        add_edge(face.a, face.b);
        add_edge(face.b, face.c);
        add_edge(face.c, face.a);
        face_removed[idx] = 1;
        face_keys.erase(CanonicalFaceKey(face.a, face.b, face.c));
    }

    std::vector<Face> new_faces;
    new_faces.reserve(boundary_edges.size());
    for (const auto &entry : boundary_edges) {
        int a = entry.first.from;
        int b = entry.first.to;
        Face face = MakeFace(a, b, new_index, vertices_, interior_point_);
        if (IsFaceDegenerate(face)) {
            continue;
        }
        FaceKey key = CanonicalFaceKey(face.a, face.b, face.c);
        if (!face_keys.insert(key).second) {
            continue;
        }
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

C3DConvexPoly C3DConvexPoly::BuildJarvisMarch(const std::vector<Vec3> &points, double epsilon) {
    TetrahedronData base = BuildInitialTetrahedron(points, epsilon);
    if (!base.valid) {
        return C3DConvexPoly();
    }

    std::vector<Vec3> hull_vertices;
    hull_vertices.reserve(points.size());
    std::vector<Face> hull_faces;
    hull_faces.reserve(points.size() * 2);
    std::unordered_map<int, int> index_map;

    auto ensure_vertex = [&](int original_index) {
        auto it = index_map.find(original_index);
        if (it != index_map.end()) {
            return it->second;
        }
        int new_index = static_cast<int>(hull_vertices.size());
        hull_vertices.push_back(points[original_index]);
        index_map.emplace(original_index, new_index);
        return new_index;
    };

    std::unordered_set<FaceKey, FaceKeyHasher> face_keys;
    face_keys.reserve(points.size() * 3);

    std::vector<DirectedEdge> edge_stack;
    edge_stack.reserve(points.size() * 6);
    std::unordered_set<DirectedEdgeKey, DirectedEdgeHasher> processed_edges;

    auto push_edge = [&](int from, int to, int previous) {
        edge_stack.push_back(DirectedEdge{from, to, previous});
    };

    auto add_face = [&](int ia, int ib, int ic, bool skip_validation) {
        int a = ensure_vertex(ia);
        int b = ensure_vertex(ib);
        int c = ensure_vertex(ic);
        Face face = BuildFace(a, b, c, hull_vertices, base.interior_point);
        if (IsFaceDegenerate(face)) {
            return;
        }
        if (!skip_validation) {
            bool valid = true;
            for (const Vec3 &point : points) {
                double distance = Dot(face.normal, point) + face.offset;
                if (distance > epsilon) {
                    valid = false;
                    break;
                }
            }
            if (!valid) {
                return;
            }
        }
        FaceKey key = CanonicalFaceKey(face.a, face.b, face.c);
        if (!face_keys.insert(key).second) {
            return;
        }
        hull_faces.push_back(face);
        push_edge(ia, ib, ic);
        push_edge(ib, ic, ia);
        push_edge(ic, ia, ib);
    };

    add_face(base.indices[0], base.indices[1], base.indices[2], true);
    add_face(base.indices[0], base.indices[3], base.indices[1], true);
    add_face(base.indices[1], base.indices[3], base.indices[2], true);
    add_face(base.indices[2], base.indices[3], base.indices[0], true);

    while (!edge_stack.empty()) {
        DirectedEdge edge = edge_stack.back();
        edge_stack.pop_back();

        DirectedEdgeKey key{edge.from, edge.to};
        if (!processed_edges.insert(key).second) {
            continue;
        }

        int candidate = FindGiftWrapCandidate(edge.from, edge.to, edge.previous, points, epsilon);
        if (candidate < 0) {
            continue;
        }

        std::size_t prev_face_count = hull_faces.size();
        add_face(edge.from, edge.to, candidate, false);
        if (hull_faces.size() == prev_face_count) {
            continue;
        }

        push_edge(edge.to, candidate, edge.from);
        push_edge(candidate, edge.from, edge.to);
    }

    std::vector<Face> filtered_faces;
    filtered_faces.reserve(hull_faces.size());
    for (const Face &face : hull_faces) {
        bool valid = true;
        for (const Vec3 &point : points) {
            double distance = Dot(face.normal, point) + face.offset;
            if (distance > epsilon) {
                valid = false;
                break;
            }
        }
        if (valid) {
            filtered_faces.push_back(face);
        }
    }

    hull_faces.swap(filtered_faces);

    return C3DConvexPoly(std::move(hull_vertices), std::move(hull_faces), base.interior_point);
}

C3DConvexPoly C3DConvexPoly::BuildGrahamScan(const std::vector<Vec3> &points, double epsilon) {
    TetrahedronData base = BuildInitialTetrahedron(points, epsilon);
    if (!base.valid) {
        return C3DConvexPoly();
    }

    C3DConvexPoly hull(std::move(base.vertices), std::move(base.faces), base.interior_point);

    std::array<int, 4> base_indices = base.indices;
    std::vector<int> order;
    order.reserve(points.size());

    auto is_base_index = [&](int idx) {
        return std::find(base_indices.begin(), base_indices.end(), idx) != base_indices.end();
    };

    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        if (!is_base_index(i)) {
            order.push_back(i);
        }
    }

    auto compute_angles = [&](const Vec3 &p) {
        Vec3 delta = p - base.interior_point;
        double radius = std::sqrt(LengthSquared(delta));
        double theta = std::atan2(delta.y, delta.x);
        double phi = std::atan2(std::sqrt(delta.x * delta.x + delta.y * delta.y), delta.z);
        return std::tuple<double, double, double>(theta, phi, radius);
    };

    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        auto lhs_angles = compute_angles(points[lhs]);
        auto rhs_angles = compute_angles(points[rhs]);
        if (std::get<0>(lhs_angles) != std::get<0>(rhs_angles)) {
            return std::get<0>(lhs_angles) < std::get<0>(rhs_angles);
        }
        if (std::get<1>(lhs_angles) != std::get<1>(rhs_angles)) {
            return std::get<1>(lhs_angles) < std::get<1>(rhs_angles);
        }
        return std::get<2>(lhs_angles) < std::get<2>(rhs_angles);
    });

    for (int idx : order) {
        hull.AddPoint(points[idx], epsilon);
    }

    return hull;
}

} // namespace c3d
