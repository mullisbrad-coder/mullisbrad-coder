#include "C3DConvexPoly.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace c3d {
namespace {

constexpr double kEpsilon = 1e-9;

double Length(const Vec3 &v) { return std::sqrt(Dot(v, v)); }
double SquaredLength(const Vec3 &v) { return Dot(v, v); }

Vec3 Normalize(const Vec3 &v) {
    double length = Length(v);
    if (length <= kEpsilon) {
        return Vec3{0.0, 0.0, 0.0};
    }
    return Vec3{v.x / length, v.y / length, v.z / length};
}

bool NearlyEqual(double a, double b, double epsilon) {
    return std::fabs(a - b) <= epsilon;
}

bool NearlyEqualVec3(const Vec3 &lhs, const Vec3 &rhs, double epsilon) {
    return std::fabs(lhs.x - rhs.x) <= epsilon && std::fabs(lhs.y - rhs.y) <= epsilon && std::fabs(lhs.z - rhs.z) <= epsilon;
}

Vec3 ComputeInteriorPoint(const std::vector<Vec3> &points) {
    Vec3 sum{0.0, 0.0, 0.0};
    if (points.empty()) {
        return sum;
    }
    for (const Vec3 &p : points) {
        sum += p;
    }
    double inv = 1.0 / static_cast<double>(points.size());
    sum *= inv;
    return sum;
}

Vec3 ComputeInteriorPoint(const std::vector<Vec3> &points, const std::vector<int> &indices) {
    std::vector<Vec3> subset;
    subset.reserve(indices.size());
    for (int idx : indices) {
        subset.push_back(points[static_cast<std::size_t>(idx)]);
    }
    return ComputeInteriorPoint(subset);
}

double SignedVolume(const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &d) {
    return Dot(Cross(b - a, c - a), d - a);
}

struct SortedFaceKey {
    int a;
    int b;
    int c;

    explicit SortedFaceKey(int ia, int ib, int ic) {
        if (ia > ib) std::swap(ia, ib);
        if (ib > ic) std::swap(ib, ic);
        if (ia > ib) std::swap(ia, ib);
        a = ia;
        b = ib;
        c = ic;
    }

    bool operator==(const SortedFaceKey &other) const {
        return a == other.a && b == other.b && c == other.c;
    }
};

struct SortedFaceKeyHasher {
    std::size_t operator()(const SortedFaceKey &key) const noexcept {
        std::size_t seed = static_cast<std::size_t>(key.a);
        seed = seed * 1315423911u + static_cast<std::size_t>(key.b);
        seed = seed * 1315423911u + static_cast<std::size_t>(key.c);
        return seed;
    }
};

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
    for (const Vec3 &existing : vertices_) {
        Vec3 diff = existing - point;
        if (SquaredLength(diff) <= epsilon * epsilon) {
            return false;
        }
    }

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

    CleanupHull(epsilon);

    return true;
}

void C3DConvexPoly::Optimize(double epsilon) {
    if (vertices_.size() >= 4) {
        std::vector<Vec3> points = vertices_;
        *this = BuildConvexPolyGraham(points, epsilon);
    } else {
        CleanupHull(epsilon);
    }
}

void C3DConvexPoly::CleanupHull(double epsilon) {
    if (faces_.empty()) {
        return;
    }

    std::vector<Face> filtered;
    filtered.reserve(faces_.size());
    std::unordered_set<SortedFaceKey, SortedFaceKeyHasher> seen;

    const double min_area = epsilon * epsilon;
    for (const Face &face : faces_) {
        const Vec3 &va = vertices_[face.a];
        const Vec3 &vb = vertices_[face.b];
        const Vec3 &vc = vertices_[face.c];
        Vec3 ab = vb - va;
        Vec3 ac = vc - va;
        Vec3 normal = Cross(ab, ac);
        double area = SquaredLength(normal);
        if (area <= min_area) {
            continue;
        }

        SortedFaceKey key(face.a, face.b, face.c);
        if (!seen.insert(key).second) {
            continue;
        }

        filtered.push_back(MakeFace(face.a, face.b, face.c, vertices_, interior_point_));
    }

    faces_.swap(filtered);

    if (faces_.empty()) {
        return;
    }

    std::vector<char> used(vertices_.size(), 0);
    for (const Face &face : faces_) {
        used[face.a] = 1;
        used[face.b] = 1;
        used[face.c] = 1;
    }

    bool all_used = true;
    for (std::size_t i = 0; i < used.size(); ++i) {
        if (!used[i]) {
            all_used = false;
            break;
        }
    }

    Vec3 new_interior = ComputeInteriorPoint(vertices_);

    if (!all_used) {
        std::vector<Vec3> new_vertices;
        new_vertices.reserve(faces_.size());
        std::vector<int> remap(vertices_.size(), -1);
        for (std::size_t i = 0; i < vertices_.size(); ++i) {
            if (used[i]) {
                remap[i] = static_cast<int>(new_vertices.size());
                new_vertices.push_back(vertices_[i]);
            }
        }

        new_interior = ComputeInteriorPoint(new_vertices);

        for (Face &face : faces_) {
            face.a = remap[face.a];
            face.b = remap[face.b];
            face.c = remap[face.c];
            face = MakeFace(face.a, face.b, face.c, new_vertices, new_interior);
        }

        vertices_.swap(new_vertices);
    }

    interior_point_ = new_interior;
}

C3DConvexPoly BuildConvexPolyJarvis(const std::vector<Vec3> &points, double epsilon) {
    if (points.size() < 4) {
        throw std::invalid_argument("Jarvis hull requires at least four non-coplanar points");
    }

    std::vector<Vec3> unique_points;
    unique_points.reserve(points.size());
    for (const Vec3 &p : points) {
        bool found = false;
        for (const Vec3 &q : unique_points) {
            if (NearlyEqualVec3(p, q, epsilon)) {
                found = true;
                break;
            }
        }
        if (!found) {
            unique_points.push_back(p);
        }
    }

    if (unique_points.size() < 4) {
        throw std::invalid_argument("Jarvis hull requires at least four unique points");
    }

    const int count = static_cast<int>(unique_points.size());

    std::array<int, 4> simplex{0, -1, -1, -1};
    for (int i = 1; i < count; ++i) {
        if (!NearlyEqualVec3(unique_points[i], unique_points[simplex[0]], epsilon)) {
            simplex[1] = i;
            break;
        }
    }
    if (simplex[1] == -1) {
        throw std::invalid_argument("Jarvis hull points are degenerate");
    }

    Vec3 base = unique_points[simplex[1]] - unique_points[simplex[0]];
    for (int i = 0; i < count; ++i) {
        if (i == simplex[0] || i == simplex[1]) continue;
        Vec3 dir = unique_points[i] - unique_points[simplex[0]];
        Vec3 cross = Cross(base, dir);
        if (Length(cross) > epsilon) {
            simplex[2] = i;
            break;
        }
    }
    if (simplex[2] == -1) {
        throw std::invalid_argument("Jarvis hull points are colinear");
    }

    Vec3 normal = Cross(unique_points[simplex[1]] - unique_points[simplex[0]], unique_points[simplex[2]] - unique_points[simplex[0]]);
    for (int i = 0; i < count; ++i) {
        if (i == simplex[0] || i == simplex[1] || i == simplex[2]) continue;
        double volume = std::fabs(SignedVolume(unique_points[simplex[0]], unique_points[simplex[1]], unique_points[simplex[2]], unique_points[i]));
        if (volume > epsilon) {
            simplex[3] = i;
            break;
        }
    }
    if (simplex[3] == -1) {
        throw std::invalid_argument("Jarvis hull points are coplanar");
    }

    Vec3 interior_point = ComputeInteriorPoint(unique_points, {simplex[0], simplex[1], simplex[2], simplex[3]});

    C3DConvexPoly poly(1.0);
    poly.vertices_.clear();
    poly.faces_.clear();

    std::vector<int> index_map(count, -1);
    for (int idx : simplex) {
        index_map[idx] = static_cast<int>(poly.vertices_.size());
        poly.vertices_.push_back(unique_points[idx]);
    }
    poly.interior_point_ = interior_point;

    auto add_face = [&](int a, int b, int c) {
        poly.faces_.push_back(poly.MakeFace(a, b, c, poly.vertices_, poly.interior_point_));
    };

    add_face(index_map[simplex[0]], index_map[simplex[1]], index_map[simplex[2]]);
    add_face(index_map[simplex[0]], index_map[simplex[2]], index_map[simplex[3]]);
    add_face(index_map[simplex[0]], index_map[simplex[3]], index_map[simplex[1]]);
    add_face(index_map[simplex[1]], index_map[simplex[3]], index_map[simplex[2]]);

    std::queue<C3DConvexPoly::EdgeKey> pending;
    auto push_edge = [&](int from, int to) {
        pending.push({from, to});
    };

    push_edge(index_map[simplex[0]], index_map[simplex[1]]);
    push_edge(index_map[simplex[1]], index_map[simplex[0]]);
    push_edge(index_map[simplex[1]], index_map[simplex[2]]);
    push_edge(index_map[simplex[2]], index_map[simplex[1]]);
    push_edge(index_map[simplex[2]], index_map[simplex[0]]);
    push_edge(index_map[simplex[0]], index_map[simplex[2]]);
    push_edge(index_map[simplex[0]], index_map[simplex[3]]);
    push_edge(index_map[simplex[3]], index_map[simplex[0]]);
    push_edge(index_map[simplex[3]], index_map[simplex[1]]);
    push_edge(index_map[simplex[1]], index_map[simplex[3]]);
    push_edge(index_map[simplex[2]], index_map[simplex[3]]);
    push_edge(index_map[simplex[3]], index_map[simplex[2]]);

    std::unordered_set<SortedFaceKey, SortedFaceKeyHasher> face_keys;
    for (const Face &face : poly.faces_) {
        face_keys.emplace(face.a, face.b, face.c);
    }

    std::unordered_map<int, int> reverse_map;
    for (int idx : simplex) {
        reverse_map[index_map[idx]] = idx;
    }

    std::unordered_set<C3DConvexPoly::EdgeKey, C3DConvexPoly::EdgeKeyHasher> processed_edges;

    while (!pending.empty()) {
        C3DConvexPoly::EdgeKey edge = pending.front();
        pending.pop();

        if (!processed_edges.insert(edge).second) {
            continue;
        }

        int a_local = edge.from;
        int b_local = edge.to;
        int a_global = reverse_map.count(a_local) ? reverse_map[a_local] : -1;
        int b_global = reverse_map.count(b_local) ? reverse_map[b_local] : -1;
        if (a_global == -1 || b_global == -1) {
            continue;
        }

        Vec3 base_vec = poly.vertices_[b_local] - poly.vertices_[a_local];
        int best = -1;
        double best_area = 0.0;

        for (int candidate = 0; candidate < count; ++candidate) {
            if (candidate == a_global || candidate == b_global) {
                continue;
            }
            Vec3 normal = Cross(base_vec, unique_points[candidate] - poly.vertices_[a_local]);
            double area = Length(normal);
            if (area <= epsilon) {
                continue;
            }

            normal = normal / area;
            double offset = -Dot(normal, poly.vertices_[a_local]);
            if (Dot(normal, poly.interior_point_) + offset > 0.0) {
                normal *= -1.0;
                offset = -offset;
            }

            bool valid = true;
            for (const Vec3 &pt : unique_points) {
                double dist = Dot(normal, pt) + offset;
                if (dist > epsilon) {
                    valid = false;
                    break;
                }
            }
            if (!valid) {
                continue;
            }

            int candidate_local = index_map[candidate];
            if (candidate_local == -1) {
                candidate_local = static_cast<int>(poly.vertices_.size());
            }
            SortedFaceKey candidate_key(a_local, b_local, candidate_local);
            if (face_keys.find(candidate_key) != face_keys.end()) {
                continue;
            }

            if (best == -1 || area > best_area + epsilon) {
                best = candidate;
                best_area = area;
            }
        }

        if (best == -1) {
            continue;
        }

        int best_local;
        if (index_map[best] == -1) {
            index_map[best] = static_cast<int>(poly.vertices_.size());
            best_local = index_map[best];
            poly.vertices_.push_back(unique_points[best]);
            reverse_map[best_local] = best;
            poly.interior_point_ = ComputeInteriorPoint(poly.vertices_);
        } else {
            best_local = index_map[best];
        }

        SortedFaceKey key(a_local, b_local, best_local);
        if (face_keys.find(key) != face_keys.end()) {
            continue;
        }

        Face face = poly.MakeFace(a_local, b_local, best_local, poly.vertices_, poly.interior_point_);
        poly.faces_.push_back(face);
        face_keys.insert(key);

        push_edge(b_local, best_local);
        push_edge(best_local, b_local);
        push_edge(best_local, a_local);
        push_edge(a_local, best_local);
    }

    poly.interior_point_ = ComputeInteriorPoint(poly.vertices_);
    poly.CleanupHull(epsilon);

    return poly;
}

C3DConvexPoly BuildConvexPolyGraham(const std::vector<Vec3> &points, double epsilon) {
    if (points.size() < 4) {
        throw std::invalid_argument("Graham hull requires at least four non-coplanar points");
    }

    std::vector<Vec3> unique_points;
    unique_points.reserve(points.size());
    for (const Vec3 &p : points) {
        bool found = false;
        for (const Vec3 &q : unique_points) {
            if (NearlyEqualVec3(p, q, epsilon)) {
                found = true;
                break;
            }
        }
        if (!found) {
            unique_points.push_back(p);
        }
    }

    if (unique_points.size() < 4) {
        throw std::invalid_argument("Graham hull requires at least four unique points");
    }

    const int count = static_cast<int>(unique_points.size());

    std::array<int, 4> simplex{0, -1, -1, -1};
    for (int i = 1; i < count; ++i) {
        if (!NearlyEqualVec3(unique_points[i], unique_points[simplex[0]], epsilon)) {
            simplex[1] = i;
            break;
        }
    }
    if (simplex[1] == -1) {
        throw std::invalid_argument("Graham hull points are degenerate");
    }

    Vec3 base = unique_points[simplex[1]] - unique_points[simplex[0]];
    for (int i = 0; i < count; ++i) {
        if (i == simplex[0] || i == simplex[1]) continue;
        Vec3 dir = unique_points[i] - unique_points[simplex[0]];
        Vec3 cross = Cross(base, dir);
        if (Length(cross) > epsilon) {
            simplex[2] = i;
            break;
        }
    }
    if (simplex[2] == -1) {
        throw std::invalid_argument("Graham hull points are colinear");
    }

    for (int i = 0; i < count; ++i) {
        if (i == simplex[0] || i == simplex[1] || i == simplex[2]) continue;
        double volume = std::fabs(SignedVolume(unique_points[simplex[0]], unique_points[simplex[1]], unique_points[simplex[2]], unique_points[i]));
        if (volume > epsilon) {
            simplex[3] = i;
            break;
        }
    }
    if (simplex[3] == -1) {
        throw std::invalid_argument("Graham hull points are coplanar");
    }

    C3DConvexPoly poly(1.0);
    poly.vertices_.clear();
    poly.faces_.clear();

    std::vector<int> remap(count, -1);
    for (int idx : simplex) {
        remap[idx] = static_cast<int>(poly.vertices_.size());
        poly.vertices_.push_back(unique_points[idx]);
    }

    poly.interior_point_ = ComputeInteriorPoint(poly.vertices_);

    auto add_face = [&](int a, int b, int c) {
        poly.faces_.push_back(poly.MakeFace(a, b, c, poly.vertices_, poly.interior_point_));
    };

    add_face(remap[simplex[0]], remap[simplex[1]], remap[simplex[2]]);
    add_face(remap[simplex[0]], remap[simplex[2]], remap[simplex[3]]);
    add_face(remap[simplex[0]], remap[simplex[3]], remap[simplex[1]]);
    add_face(remap[simplex[1]], remap[simplex[3]], remap[simplex[2]]);

    struct SphericalKey {
        double azimuth;
        double elevation;
        double radius;
        int index;
    };

    Vec3 centroid = ComputeInteriorPoint(unique_points);
    std::vector<SphericalKey> order;
    order.reserve(unique_points.size());

    for (int i = 0; i < count; ++i) {
        if (remap[i] != -1) {
            continue;
        }
        Vec3 diff = unique_points[i] - centroid;
        double azimuth = std::atan2(diff.y, diff.x);
        double radius_xy = std::sqrt(diff.x * diff.x + diff.y * diff.y);
        double elevation = std::atan2(diff.z, radius_xy);
        double radius = Length(diff);
        order.push_back({azimuth, elevation, radius, i});
    }

    std::sort(order.begin(), order.end(), [](const SphericalKey &lhs, const SphericalKey &rhs) {
        if (!NearlyEqual(lhs.azimuth, rhs.azimuth, 1e-12)) {
            return lhs.azimuth < rhs.azimuth;
        }
        if (!NearlyEqual(lhs.elevation, rhs.elevation, 1e-12)) {
            return lhs.elevation < rhs.elevation;
        }
        return lhs.radius < rhs.radius;
    });

    for (const SphericalKey &key : order) {
        poly.AddPoint(unique_points[key.index], epsilon);
    }

    poly.CleanupHull(epsilon);

    return poly;
}

} // namespace c3d

