#include "C3DConvexPoly.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace c3d {
namespace {

constexpr double kEpsilon = 1e-9;

double Length(const Vec3 &v) {
    return std::sqrt(Dot(v, v));
}

Vec3 Normalize(const Vec3 &v) {
    double length = Length(v);
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

C3DConvexPoly::C3DConvexPoly(double half_extent, HorizonTriangulation method) : method_(method) {
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

Face C3DConvexPoly::MakeFace(int a, int b, int c, const std::vector<Vec3> &vertices, const Vec3 &interior_point) {
    return BuildFace(a, b, c, vertices, interior_point);
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

    std::unordered_set<EdgeKey, EdgeKeyHasher> boundary_edges;
    auto add_edge = [&](int from, int to) {
        EdgeKey key{from, to};
        EdgeKey opposite{to, from};
        if (boundary_edges.erase(opposite) == 0) {
            boundary_edges.insert(key);
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

    std::vector<EdgeKey> edge_list;
    edge_list.reserve(boundary_edges.size());
    for (const auto &edge : boundary_edges) {
        edge_list.push_back(edge);
    }

    std::vector<Face> new_faces = BuildHorizonFaces(new_index, edge_list);
    if (new_faces.empty()) {
        new_faces.reserve(edge_list.size());
        for (const EdgeKey &edge : edge_list) {
            new_faces.push_back(MakeFace(edge.from, edge.to, new_index, vertices_, interior_point_));
        }
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

std::vector<C3DConvexPoly::ProjectedVertex> C3DConvexPoly::ProjectComponent(
    const std::vector<int> &component_vertices,
    const std::vector<int> &component_edges,
    int new_index,
    const std::vector<EdgeKey> &edges) const {
    std::vector<ProjectedVertex> projected;
    if (component_vertices.empty()) {
        return projected;
    }

    const Vec3 &origin = vertices_[new_index];
    Vec3 normal{0.0, 0.0, 0.0};
    for (int edge_idx : component_edges) {
        const EdgeKey &edge = edges[edge_idx];
        Vec3 a = vertices_[edge.from] - origin;
        Vec3 b = vertices_[edge.to] - origin;
        normal += Cross(a, b);
    }

    if (Length(normal) <= kEpsilon && component_vertices.size() >= 2) {
        Vec3 a = vertices_[component_vertices[0]] - origin;
        Vec3 b = vertices_[component_vertices[1]] - origin;
        normal = Cross(a, b);
    }
    if (Length(normal) <= kEpsilon) {
        normal = Vec3{0.0, 0.0, 1.0};
    }
    normal = Normalize(normal);

    Vec3 u{0.0, 0.0, 0.0};
    for (int idx : component_vertices) {
        Vec3 candidate = vertices_[idx] - origin;
        candidate -= normal * Dot(candidate, normal);
        if (Length(candidate) > kEpsilon) {
            u = Normalize(candidate);
            break;
        }
    }
    if (Length(u) <= kEpsilon) {
        Vec3 arbitrary = std::fabs(normal.x) < 0.9 ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
        Vec3 projected_arbitrary = arbitrary - normal * Dot(arbitrary, normal);
        if (Length(projected_arbitrary) <= kEpsilon) {
            arbitrary = Vec3{0.0, 0.0, 1.0};
            projected_arbitrary = arbitrary - normal * Dot(arbitrary, normal);
        }
        u = Normalize(projected_arbitrary);
    }
    Vec3 v = Normalize(Cross(normal, u));
    if (Length(v) <= kEpsilon) {
        Vec3 arbitrary = Vec3{0.0, 0.0, 1.0};
        Vec3 projected_arbitrary = arbitrary - normal * Dot(arbitrary, normal);
        u = Normalize(projected_arbitrary);
        v = Normalize(Cross(normal, u));
    }

    projected.reserve(component_vertices.size());
    for (int idx : component_vertices) {
        Vec3 relative = vertices_[idx] - origin;
        double x = Dot(relative, u);
        double y = Dot(relative, v);
        projected.push_back(ProjectedVertex{idx, x, y});
    }
    return projected;
}

double C3DConvexPoly::Orientation2D(const ProjectedVertex &a,
                                    const ProjectedVertex &b,
                                    const ProjectedVertex &c) {
    double ux = b.x - a.x;
    double uy = b.y - a.y;
    double vx = c.x - a.x;
    double vy = c.y - a.y;
    return ux * vy - uy * vx;
}

double C3DConvexPoly::DistanceSq2D(const ProjectedVertex &a, const ProjectedVertex &b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return dx * dx + dy * dy;
}

std::vector<int> C3DConvexPoly::JarvisMarchHull(const std::vector<ProjectedVertex> &points) const {
    std::vector<int> hull;
    const int n = static_cast<int>(points.size());
    if (n == 0) {
        return hull;
    }
    if (n <= 2) {
        for (const auto &p : points) {
            hull.push_back(p.index);
        }
        return hull;
    }

    int leftmost = 0;
    for (int i = 1; i < n; ++i) {
        if (points[i].x < points[leftmost].x - kEpsilon ||
            (std::fabs(points[i].x - points[leftmost].x) <= kEpsilon && points[i].y < points[leftmost].y)) {
            leftmost = i;
        }
    }

    int p = leftmost;
    do {
        hull.push_back(points[p].index);
        int q = (p + 1) % n;
        if (q == p) {
            q = (p + 2) % n;
        }
        for (int r = 0; r < n; ++r) {
            if (r == p || r == q) {
                continue;
            }
            double orient = Orientation2D(points[p], points[q], points[r]);
            if (orient > kEpsilon ||
                (std::fabs(orient) <= kEpsilon &&
                 DistanceSq2D(points[p], points[r]) > DistanceSq2D(points[p], points[q]))) {
                q = r;
            }
        }
        p = q;
    } while (p != leftmost);

    return hull;
}

std::vector<int> C3DConvexPoly::GrahamScanHull(const std::vector<ProjectedVertex> &points) const {
    std::vector<int> hull_indices;
    if (points.empty()) {
        return hull_indices;
    }
    if (points.size() <= 2) {
        for (const auto &p : points) {
            hull_indices.push_back(p.index);
        }
        return hull_indices;
    }

    std::vector<ProjectedVertex> pts = points;
    std::sort(pts.begin(), pts.end(), [](const ProjectedVertex &a, const ProjectedVertex &b) {
        if (std::fabs(a.x - b.x) > kEpsilon) {
            return a.x < b.x;
        }
        if (std::fabs(a.y - b.y) > kEpsilon) {
            return a.y < b.y;
        }
        return a.index < b.index;
    });

    pts.erase(std::unique(pts.begin(), pts.end(), [](const ProjectedVertex &a, const ProjectedVertex &b) {
                  return std::fabs(a.x - b.x) <= kEpsilon && std::fabs(a.y - b.y) <= kEpsilon;
              }),
              pts.end());

    if (pts.size() <= 2) {
        for (const auto &p : pts) {
            hull_indices.push_back(p.index);
        }
        return hull_indices;
    }

    std::vector<ProjectedVertex> lower;
    for (const auto &pt : pts) {
        while (lower.size() >= 2) {
            double orient = Orientation2D(lower[lower.size() - 2], lower.back(), pt);
            if (orient > kEpsilon) {
                break;
            }
            lower.pop_back();
        }
        lower.push_back(pt);
    }

    std::vector<ProjectedVertex> upper;
    for (auto it = pts.rbegin(); it != pts.rend(); ++it) {
        while (upper.size() >= 2) {
            double orient = Orientation2D(upper[upper.size() - 2], upper.back(), *it);
            if (orient > kEpsilon) {
                break;
            }
            upper.pop_back();
        }
        upper.push_back(*it);
    }

    lower.pop_back();
    upper.pop_back();

    std::vector<ProjectedVertex> hull;
    hull.reserve(lower.size() + upper.size());
    hull.insert(hull.end(), lower.begin(), lower.end());
    hull.insert(hull.end(), upper.begin(), upper.end());

    for (const auto &pt : hull) {
        hull_indices.push_back(pt.index);
    }
    return hull_indices;
}

std::vector<Face> C3DConvexPoly::BuildHorizonFaces(int new_index, const std::vector<EdgeKey> &edges) const {
    std::vector<Face> faces;
    if (edges.empty()) {
        return faces;
    }

    if (method_ == HorizonTriangulation::kEdgeFan) {
        faces.reserve(edges.size());
        for (const EdgeKey &edge : edges) {
            faces.push_back(MakeFace(edge.from, edge.to, new_index, vertices_, interior_point_));
        }
        return faces;
    }

    std::unordered_map<int, std::vector<int>> vertex_to_edges;
    vertex_to_edges.reserve(edges.size() * 2);
    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
        const EdgeKey &edge = edges[i];
        vertex_to_edges[edge.from].push_back(i);
        vertex_to_edges[edge.to].push_back(i);
    }

    std::unordered_set<int> visited_vertices;
    std::vector<char> edge_used(edges.size(), 0);

    for (const auto &entry : vertex_to_edges) {
        int start_vertex = entry.first;
        if (!visited_vertices.insert(start_vertex).second) {
            continue;
        }

        std::queue<int> queue;
        std::vector<int> component_vertices;
        std::vector<int> component_edges;

        queue.push(start_vertex);
        component_vertices.push_back(start_vertex);

        while (!queue.empty()) {
            int v = queue.front();
            queue.pop();

            const auto &adj_edges = vertex_to_edges.at(v);
            for (int edge_idx : adj_edges) {
                if (!edge_used[edge_idx]) {
                    edge_used[edge_idx] = 1;
                    component_edges.push_back(edge_idx);
                }

                const EdgeKey &edge = edges[edge_idx];
                int neighbor = edge.from == v ? edge.to : edge.from;
                if (visited_vertices.insert(neighbor).second) {
                    component_vertices.push_back(neighbor);
                    queue.push(neighbor);
                }
            }
        }

        if (component_vertices.size() < 3) {
            continue;
        }

        std::vector<ProjectedVertex> projected =
            ProjectComponent(component_vertices, component_edges, new_index, edges);
        if (projected.empty()) {
            continue;
        }

        std::vector<int> horizon_order;
        if (method_ == HorizonTriangulation::kJarvisMarch) {
            horizon_order = JarvisMarchHull(projected);
        } else {
            horizon_order = GrahamScanHull(projected);
        }

        if (horizon_order.size() < 3) {
            for (int edge_idx : component_edges) {
                const EdgeKey &edge = edges[edge_idx];
                faces.push_back(MakeFace(edge.from, edge.to, new_index, vertices_, interior_point_));
            }
            continue;
        }

        for (size_t i = 0; i < horizon_order.size(); ++i) {
            int a = horizon_order[i];
            int b = horizon_order[(i + 1) % horizon_order.size()];
            faces.push_back(MakeFace(a, b, new_index, vertices_, interior_point_));
        }
    }

    return faces;
}

} // namespace c3d
