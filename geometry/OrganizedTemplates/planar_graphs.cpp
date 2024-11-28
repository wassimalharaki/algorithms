/**
 * 22:47:20 11/23/24
 * planar_graphs
 */
// ./algorithms/geometry/OrganizedTemplates/planar_graphs.cpp
#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;

/**
 * facts:
 *
 * n: vertices.
 * m: edges.
 * f: faces.
 * n - m + f = 2
 *
 * More generally:
 * n - m + f = 1 + k, where k is the number of connected components.
 *
 * if n >= 3, the max number of edges in a planar graph with n vertices is 3 * n - 6.
 * as a direct consequence, the max number of faces is 2 * n - 4.
 *
 * Every planar graph has a vertex of degree 5 or less.
 */

// planar graph from points
namespace pnt {
    struct Point {
        int64_t x, y;

        Point(int64_t x_, int64_t y_): x(x_), y(y_) {}

        Point operator - (const Point & p) const {
            return Point(x - p.x, y - p.y);
        }

        int64_t cross (const Point & p) const {
            return x * p.y - y * p.x;
        }

        int64_t cross (const Point & p, const Point & q) const {
            return (p - *this).cross(q - *this);
        }

        int half () const {
            return int(y < 0 || (y == 0 && x < 0));
        }
    };

    std::vector<std::vector<size_t>> find_faces(std::vector<Point> vertices, std::vector<std::vector<size_t>> adj) {
        size_t n = vertices.size();
        std::vector<std::vector<char>> used(n);
        for (size_t i = 0; i < n; i++) {
            used[i].resize(adj[i].size());
            used[i].assign(adj[i].size(), 0);
            auto compare = [&](size_t l, size_t r) {
                Point pl = vertices[l] - vertices[i];
                Point pr = vertices[r] - vertices[i];
                if (pl.half() != pr.half())
                    return pl.half() < pr.half();
                return pl.cross(pr) > 0;
            };
            std::sort(adj[i].begin(), adj[i].end(), compare);
        }
        std::vector<std::vector<size_t>> faces;
        for (size_t i = 0; i < n; i++) {
            for (size_t edge_id = 0; edge_id < adj[i].size(); edge_id++) {
                if (used[i][edge_id]) {
                    continue;
                }
                std::vector<size_t> face;
                size_t v = i;
                size_t e = edge_id;
                while (!used[v][e]) {
                    used[v][e] = true;
                    face.push_back(v);
                    size_t u = adj[v][e];
                    size_t e1 = std::lower_bound(adj[u].begin(), adj[u].end(), v, [&](size_t l, size_t r) {
                        Point pl = vertices[l] - vertices[u];
                        Point pr = vertices[r] - vertices[u];
                        if (pl.half() != pr.half())
                            return pl.half() < pr.half();
                        return pl.cross(pr) > 0;
                    }) - adj[u].begin() + 1;
                    if (e1 == adj[u].size()) {
                        e1 = 0;
                    }
                    v = u;
                    e = e1;
                }
                std::reverse(face.begin(), face.end());
                int sign = 0;
                for (size_t j = 0; j < face.size(); j++) {
                    size_t j1 = (j + 1) % face.size();
                    size_t j2 = (j + 2) % face.size();
                    int64_t val = vertices[face[j]].cross(vertices[face[j1]], vertices[face[j2]]);
                    if (val > 0) {
                        sign = 1;
                        break;
                    } else if (val < 0) {
                        sign = -1;
                        break;
                    }
                }
                if (sign <= 0) {
                    faces.insert(faces.begin(), face);
                } else {
                    faces.emplace_back(face);
                }
            }
        }
        return faces;
    }
}

// planar graph from segments
namespace seg {
    using dbl = long double;

    const dbl eps = 1e-9;

    struct Point {
        dbl x, y;

        Point(){}
        Point(dbl x_, dbl y_): x(x_), y(y_) {}

        Point operator * (dbl d) const {
            return Point(x * d, y * d);
        }

        Point operator + (const Point & p) const {
            return Point(x + p.x, y + p.y);
        }

        Point operator - (const Point & p) const {
            return Point(x - p.x, y - p.y);
        }

        dbl cross (const Point & p) const {
            return x * p.y - y * p.x;
        }

        dbl cross (const Point & p, const Point & q) const {
            return (p - *this).cross(q - *this);
        }

        dbl dot (const Point & p) const {
            return x * p.x + y * p.y;
        }

        dbl dot (const Point & p, const Point & q) const {
            return (p - *this).dot(q - *this);
        }

        bool operator < (const Point & p) const {
            if (fabs(x - p.x) < eps) {
                if (fabs(y - p.y) < eps) {
                    return false;
                } else {
                    return y < p.y;
                }
            } else {
                return x < p.x;
            }
        }

        bool operator == (const Point & p) const {
            return fabs(x - p.x) < eps && fabs(y - p.y) < eps;
        }

        bool operator >= (const Point & p) const {
            return !(*this < p);
        }
    };

    struct Line{
        Point p[2];

        Line(Point l, Point r){p[0] = l; p[1] = r;}
        Point& operator [](const int & i){return p[i];}
        const Point& operator[](const int & i)const{return p[i];}
        Line(const Line & l){
            p[0] = l.p[0]; p[1] = l.p[1];
        }
        Point getOrth()const{
            return Point(p[1].y - p[0].y, p[0].x - p[1].x);
        }
        bool hasPointLine(const Point & t)const{
            return std::fabs(p[0].cross(p[1], t)) < eps;
        }
        bool hasPointSeg(const Point & t)const{
            return hasPointLine(t) && t.dot(p[0], p[1]) < eps;
        }
    };

    std::vector<Point> interLineLine(Line l1, Line l2){
        if(std::fabs(l1.getOrth().cross(l2.getOrth())) < eps){
            if(l1.hasPointLine(l2[0]))return {l1[0], l1[1]};
            else return {};
        }
        Point u = l2[1] - l2[0];
        Point v = l1[1] - l1[0];
        dbl s = u.cross(l2[0] - l1[0])/u.cross(v);
        return {Point(l1[0] + v * s)};
    }

    std::vector<Point> interSegSeg(Line l1, Line l2){
        if (l1[0] == l1[1]) {
            if (l2[0] == l2[1]) {
                if (l1[0] == l2[0])
                    return {l1[0]};
                else
                    return {};
            } else {
                if (l2.hasPointSeg(l1[0]))
                    return {l1[0]};
                else
                    return {};
            }
        }
        if (l2[0] == l2[1]) {
            if (l1.hasPointSeg(l2[0]))
                return {l2[0]};
            else
                return {};
        }
        auto li = interLineLine(l1, l2);
        if (li.empty())
            return li;
        if (li.size() == 2) {
            if (l1[0] >= l1[1])
                std::swap(l1[0], l1[1]);
            if (l2[0] >= l2[1])
                std::swap(l2[0], l2[1]);
            std::vector<Point> res(2);
            if (l1[0] < l2[0])
                res[0] = l2[0];
            else
                res[0] = l1[0];
            if (l1[1] < l2[1])
                res[1] = l1[1];
            else
                res[1] = l2[1];
            if (res[0] == res[1])
                res.pop_back();
            if (res.size() == 2u && res[1] < res[0])
                return {};
            else
                return res;
        }
        Point cand = li[0];
        if (l1.hasPointSeg(cand) && l2.hasPointSeg(cand))
            return {cand};
        else
            return {};
    }

    std::pair<std::vector<Point>, std::vector<std::vector<size_t>>> build_graph(std::vector<Line> segments) {
        std::vector<Point> p;
        std::vector<std::vector<size_t>> adj;
        std::map<std::pair<int64_t, int64_t>, size_t> point_id;
        auto get_point_id = [&](Point pt) {
            auto repr = std::make_pair(
                    int64_t(std::round(pt.x * 1000000000) + 1e-6),
                    int64_t(std::round(pt.y * 1000000000) + 1e-6)
            );
            if (!point_id.count(repr)) {
                adj.emplace_back();
                size_t id = point_id.size();
                point_id[repr] = id;
                p.push_back(pt);
                return id;
            } else {
                return point_id[repr];
            }
        };
        for (size_t i = 0; i < segments.size(); i++) {
            std::vector<size_t> curr = {
                    get_point_id(segments[i][0]),
                    get_point_id(segments[i][1])
            };
            for (size_t j = 0; j < segments.size(); j++) {
                if (i == j)
                    continue;
                auto inter = interSegSeg(segments[i], segments[j]);
                for (auto pt: inter) {
                    curr.push_back(get_point_id(pt));
                }
            }
            std::sort(curr.begin(), curr.end(), [&](size_t l, size_t r) { return p[l] < p[r]; });
            curr.erase(std::unique(curr.begin(), curr.end()), curr.end());
            for (size_t j = 0; j + 1 < curr.size(); j++) {
                adj[curr[j]].push_back(curr[j + 1]);
                adj[curr[j + 1]].push_back(curr[j]);
            }
        }
        for (size_t i = 0; i < adj.size(); i++) {
            std::sort(adj[i].begin(), adj[i].end());
            // removing edges that were added multiple times
            adj[i].erase(std::unique(adj[i].begin(), adj[i].end()), adj[i].end());
        }
        return {p, adj};
    }
}
