/**
 * 16:21:39 9/3/24
 * shortest_path
 */
// ./ICPC/Geometry/Templates/shortest_path.cpp
// https://github.com/kevinzg/spholes/blob/master/src/MainWindow.h
#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;
#define int long long
#define uint unsigned int
#define double long double
#define fast ios_base::sync_with_stdio(false);cin.tie(NULL);
#define nl '\n'
#define all(v) v.begin(), v.end()
#define NO (cout << "NO" << nl);
#define YES (cout << "YES" << nl);
#define F first
#define S second
#define MOD 1000000007ll
#define pii pair<int, int>
#define P complex<int>
#define X real()
#define Y imag()
#define vec vector
#define LT int T; cin >> T; while (T--)
#define ftype double
template<typename T>
using vec2d = vector<vector<T>>;
template<typename T>
using ordered_set = tree<T, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;
using indexed_set = tree<int, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;

const ftype EPS = 1e-9, INF = 1e9, PI = M_PI;

bool approx(ftype a, ftype b) {
    return abs(a - b) < EPS;
}
bool lessThan(ftype a, ftype b) {
    return a - b < EPS;
}

struct PolarPoint;

struct Point: public complex<ftype> {
    Point() = default;
    Point(ftype x, ftype y): complex<ftype>(x, y) {}
    Point(const complex<ftype>& c): Point(c.X, c.Y) {}
    inline ftype x() const { return this->X; }
    inline ftype y() const { return this->Y; }
    inline PolarPoint toPolarPoint() const;

    inline bool operator<(const Point& a) const {
        if (approx(x(), a.x())) return lessThan(y(), a.y());
        return lessThan(x(), a.x());
    }

    inline ftype cross(const Point& b) const {
        return x() * b.y() - y() * b.x();
    }

    inline ftype dot(const Point& b) const {
        return x() * b.x() + y() * b.y();
    }

};
istream& operator>>(istream& in, P& p) {
    ftype x, y;
    in >> x >> y;
    p = {x, y};
    return in;
}

struct PolarPoint: public complex<ftype> {
    PolarPoint() = default;
    PolarPoint(ftype angle, ftype radius): complex<ftype>(angle, radius) {}
    inline ftype angle() const { return this->X; }
    inline ftype radius() const { return this->Y; }
    inline Point toPoint() const {
        return {cos(angle()) * radius(), sin(angle()) * radius()};
    }
};

PolarPoint Point::toPolarPoint() const {
    ftype angle = atan2(this->y(), this->x());
    return {angle < 0.0 ? angle + 2 * PI : angle, hypot(this->x(), this->y())};
}

struct Path: public vector<Point> {
};

struct LineSegment: public pair<Point, Point> {
    LineSegment() = default;
    LineSegment(const Point& a, const Point& b): pair<Point, Point>(a, b) {}

    enum IntersectionMode {
        NoIntersection,
        PointIntersection,
        SegmentIntersection
    };

    IntersectionMode intersection(const LineSegment& b, Point& p, Point& q) {
        Point dirA = S - F;
        Point dirB = b.S - b.F;
        Point dirAB = b.F - F;

        ftype crossAB = dirA.cross(dirB);
        ftype crossABB = dirAB.cross(dirB);
        if (approx(crossAB, 0) and approx(crossABB, 0)) {
            Point x = dirA / dirA.dot(dirA);
            ftype s = dirAB.dot(x);
            ftype t = s + dirB.dot(x);

            if (s > t) ::swap(s, t);

            s = max(s, 0.0L);
            t = min(t, 1.0L);

            if (s > 1.0 + EPS or t < 0.0 - EPS) return NoIntersection;

            p = F + dirA * s;
            q = F + dirA * t;

            return SegmentIntersection;
        }

        if (approx(crossAB, 0.0)) {
            return NoIntersection;
        }

        ftype s = dirAB.cross(dirB) / crossAB;
        ftype t = dirAB.cross(dirA) / crossAB;

        if (s >= 0.0 - EPS && s <= 1.0 + EPS and t >= 0.0 - EPS and t <= 1.0 + EPS) {
            p = F + dirA * s;
            q = b.F + dirB * t;

            return PointIntersection;
        }

        return NoIntersection;
    }
};

namespace std {
    template<>
    struct hash<Point> {
        size_t operator()(const Point& point) const {
            size_t seed = 0;
            seed ^= hash<ftype>()(point.x()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= hash<ftype>()(point.y()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };
}

struct Polygon: public vector<Point> {
    Polygon() = default;
    Polygon(int n): vec<Point>(n) {}
};

template<class V>
struct Graph {
    map<V, vec<V>> adjList;
    Graph() = default;
    void addEdge(const V& a, const V& b) {
        adjList[a].push_back(b);
    }
    inline auto begin() const {
        return adjList.begin();
    }
    inline auto end() const {
        return adjList.end();
    }
    inline const auto& getEdges(const V& p) const {
        return adjList.at(p);
    }
};

struct VisibilityGraph {
    struct PointRef;

    static vec<Point> visibleVertices(const Point &pivot, const vec<Polygon> &obstacles, PointRef pivotRef) {
        // Vector to store the points.
        std::vector<std::pair<PointRef, PolarPoint>> points;

        // Counter of all the vertices in the obstacles.
        size_t globalPointCounter = 0;

        // Iterate on all the obstacles and their vertices to add them to the _points_ vector.
        for (auto it = obstacles.begin(); it != obstacles.end(); ++it)
        {
            const Polygon &polygon = *it;
            for (auto pointIt = polygon.begin(); pointIt != polygon.end(); ++pointIt)
            {
                // We add the points in polar coordinates, translating them so _pivot_ becomes the origin.
                PolarPoint polarPoint = Point(*pointIt - pivot).toPolarPoint();
                points.push_back(std::make_pair(PointRef({static_cast<size_t>(std::distance(obstacles.begin(), it)),
                                                          static_cast<size_t>(std::distance(polygon.begin(), pointIt)),
                                                          globalPointCounter++}),
                                                polarPoint));
            }
        }

        // If the pivot is an extern point (i.e. not part of the obstacles) update its reference so calculations are easier later on.
        if (pivotRef.polygonId == obstacles.size())
            pivotRef.globalPointId = pivotRef.localPointId = globalPointCounter;

        // Sort the points by their angle, then radius.
        std::sort(points.begin(), points.end(), [](const auto &a, const auto &b) {
            if (a.second.angle() != b.second.angle())
                return a.second.angle() < b.second.angle();

            if (a.second.radius() != b.second.radius())
                return a.second.radius() < b.second.radius();

            return a.first.globalPointId < b.first.globalPointId;
        });

        // This is an index to find a point on the sorted vector by its global point id.
        std::vector<size_t> sortedPointsIndices(points.size());
        for (auto pointIt = points.begin(); pointIt != points.end(); ++pointIt)
            sortedPointsIndices[pointIt->first.globalPointId] = static_cast<size_t>(std::distance(points.begin(), pointIt));

        // Vector to store references to line segments.
        std::vector<LineSegmentRef> lineSegments(points.size());

        // Build the _lineSegments_ vector iterating over all points.
        for (auto pointIt = points.begin(); pointIt != points.end(); ++pointIt)
        {
            size_t polarPointId = static_cast<size_t>(std::distance(points.begin(), pointIt));

            size_t cwLineSegment = pointIt->first.globalPointId;
            size_t ccwLineSegment = pointIt->first.localPointId > 0 ?
                                    pointIt->first.globalPointId - 1 :
                                    pointIt->first.globalPointId + obstacles[pointIt->first.polygonId].size() - 1;

            lineSegments[cwLineSegment].polarPointAId = polarPointId;
            lineSegments[ccwLineSegment].polarPointBId = polarPointId;

            lineSegments[cwLineSegment].lineSegmentId = cwLineSegment;
        }

        // From now on the pivot will be the origin.
        const Point origin(0, 0);
        LineSegment sweepLine(origin, Point(INF, 0));

        // Line segment comparator based on the distance from the pivot to each line segment throug the sweep line.
        auto lineSegmentTreeComp = [&points, &sweepLine, &origin](const LineSegmentRef &a, const LineSegmentRef &b)
        {
            LineSegment segmentA(points[a.polarPointAId].second.toPoint(), points[a.polarPointBId].second.toPoint());
            LineSegment segmentB(points[b.polarPointAId].second.toPoint(), points[b.polarPointBId].second.toPoint());

            Point p, q;

            ftype distanceA = INF;
            ftype distanceB = INF;

            if (sweepLine.intersection(segmentA, p, q) == LineSegment::PointIntersection)
                distanceA = std::hypot(p.x(), p.y());
            if (sweepLine.intersection(segmentB, p, q) == LineSegment::PointIntersection)
                distanceB = std::hypot(p.x(), p.y());

            if (!approx(distanceA, distanceB))
                return distanceA < distanceB;

            // If the distances are the same (e.g. the line segments share a vertex on the sweep line)
            // advance the sweep line a little bit and compare again.
            ftype nextAngle = sweepLine.second.toPolarPoint().angle() + 0.01;
            LineSegment nextSweepLine(origin, PolarPoint(nextAngle, INF).toPoint());

            distanceA = INF;
            distanceB = INF;

            if (nextSweepLine.intersection(segmentA, p, q) == LineSegment::PointIntersection)
                distanceA = std::hypot(p.x(), p.y());
            if (nextSweepLine.intersection(segmentB, p, q) == LineSegment::PointIntersection)
                distanceB = std::hypot(p.x(), p.y());

            if (!approx(distanceA, distanceB))
                return distanceA < distanceB;

            return a.lineSegmentId < b.lineSegmentId;
        };

        // Line segment tree to store the line segments using the previous comparator.
        typedef std::set<LineSegmentRef, decltype(lineSegmentTreeComp)> LineSegmentTree;
        LineSegmentTree lineSegmentTree(lineSegmentTreeComp);
        std::vector<typename LineSegmentTree::iterator> lineSegmentIterators(lineSegments.size(), lineSegmentTree.end());

        // Add line segments that intersect the current sweep line.
        for (auto lineSegmentIt = lineSegments.begin(); lineSegmentIt != lineSegments.end(); ++lineSegmentIt)
        {
            const Point pointA = points[lineSegmentIt->polarPointAId].second.toPoint();
            const Point pointB = points[lineSegmentIt->polarPointBId].second.toPoint();

            Point p, dummyPointQ;

            LineSegment::IntersectionMode intersection =
                    sweepLine.intersection(LineSegment(pointA, pointB), p, dummyPointQ);

            if (intersection != LineSegment::IntersectionMode::PointIntersection)
                continue;

            ftype cpA = pointA.cross(p);
            ftype cpB = pointB.cross(p);

            // This means the line segment properly intersects the sweep line, so it has to be added to the tree.
            if (cpA * cpB < 0)
            {
                size_t lineSegmentId = static_cast<size_t>(std::distance(lineSegments.begin(), lineSegmentIt));
                lineSegmentIterators[lineSegmentId] = lineSegmentTree.insert(*lineSegmentIt).first;
            }
        }

        // Vector to store the visible points, it's the vector to be returned at the end.
        std::vector<Point> visiblePoints;

        // Iterate on all points to check if they are visible.
        for (auto pointIt = points.begin(); pointIt != points.end(); ++pointIt)
        {
            const PointRef &pointRef = pointIt->first;
            const PolarPoint&point = pointIt->second;

            // Skip pivot.
            if (pointRef.globalPointId == pivotRef.globalPointId)
                continue;

            // Move _sweepLine_ to be the ray from the origin to _point_.
            sweepLine.second = PolarPoint(point.angle(), INF).toPoint();

            // This determines if the point is visible or not.
            auto visible = [&]()
            {
                // Do a special check if _point_ and the pivot are vertices of the same polygon.
                if (pointRef.polygonId == pivotRef.polygonId)
                {
                    // Next and previous points to the pivot. Clockwise.
                    size_t pivotNextPoint = pivotRef.localPointId == obstacles[pivotRef.polygonId].size() - 1 ?
                                            pivotRef.globalPointId - (obstacles[pivotRef.polygonId].size() - 1) : pivotRef.globalPointId + 1;
                    size_t pivotPrevPoint = pivotRef.localPointId == 0 ?
                                            pivotRef.globalPointId + (obstacles[pivotRef.polygonId].size() - 1) : pivotRef.globalPointId - 1;

                    // If the two vertices are adyacent then the point is visble.
                    if (pointRef.globalPointId == pivotNextPoint || pointRef.globalPointId == pivotPrevPoint)
                        return true;

                    // If the point is between the two incident edges to the pivot then the point is not visible.
                    ftype prevAngle = points[sortedPointsIndices[pivotPrevPoint]].second.angle();
                    ftype nextAngle = points[sortedPointsIndices[pivotNextPoint]].second.angle();
                    ftype currAngle = point.angle();

                    nextAngle -= prevAngle;
                    currAngle -= prevAngle;

                    nextAngle = nextAngle < 0.0 ? nextAngle + 2 * PI : nextAngle;
                    currAngle = currAngle < 0.0 ? currAngle + 2 * PI : currAngle;

                    if (nextAngle > currAngle)
                        return false;
                }

                // If the point before this one shares the same angle with this (i.e. both are on the sweep line) then return false.
                // Note that the point might be visible, but for the shortest path problem this isn't really important since the previous point
                // will be connected to this one.
                // The right thing to do here would be to check if there is an edge on the line segments tree between
                // the previous point radius and this point radius.
                if (pointIt != points.begin() && pointIt->second.angle() == (pointIt - 1)->second.angle())
                    return false;

                // If there are no edges on the line segments tree then the point is visible
                if (lineSegmentTree.empty())
                    return true;

                // These are the points of the first line segment in the line segment tree.
                auto pointAIt = points.begin() + lineSegmentTree.begin()->polarPointAId;
                auto pointBIt = points.begin() + lineSegmentTree.begin()->polarPointBId;

                // If the point is one of the front edge ones then mark it as visible.
                if (pointAIt == pointIt || pointBIt == pointIt)
                    return true;

                // Just check if the front edge on the tree intersects the vision line from the origin to the point.
                LineSegment frontEdge(pointAIt->second.toPoint(), pointBIt->second.toPoint());
                LineSegment visionLine(origin, point.toPoint());

                Point p, dummyPointQ;

                // If the point of intersection is _point_ then it is visible.
                if (frontEdge.intersection(visionLine, p, dummyPointQ) == LineSegment::PointIntersection)
                    return approx(p.x(), point.toPoint().x()) && approx(p.y(), point.toPoint().y());

                return true;
            };

            // Check if the point is visible.
            if (visible())
                visiblePoints.push_back(obstacles[pointRef.polygonId][pointRef.localPointId]);

            // Find the line segments incident to _point_.
            size_t cwLineSegmentId = pointRef.globalPointId;
            size_t ccwLineSegmentId = pointRef.localPointId > 0 ?
                                      pointRef.globalPointId - 1 :
                                      pointRef.globalPointId + obstacles[pointRef.polygonId].size() - 1;

            const LineSegmentRef &cwLineSegment = lineSegments[cwLineSegmentId];
            const LineSegmentRef &ccwLineSegment = lineSegments[ccwLineSegmentId];

            // This adds or removes the line segment to the line segment tree.
            auto processLineSegment = [&](const LineSegmentRef& lineSegment) mutable
            {
                // If the line segment is incident to the pivot then do nothing.
                if (points[lineSegment.polarPointAId].first.globalPointId == pivotRef.globalPointId ||
                    points[lineSegment.polarPointBId].first.globalPointId == pivotRef.globalPointId)
                    return;

                auto &lineSegmentIt = lineSegmentIterators[lineSegment.lineSegmentId];

                Point otherPoint = points[lineSegment.polarPointAId].first.globalPointId == pointRef.globalPointId ?
                                   points[lineSegment.polarPointBId].second.toPoint() : points[lineSegment.polarPointAId].second.toPoint();

                ftype cp = otherPoint.cross(point.toPoint());

                if (cp < 0 && lineSegmentIt == lineSegmentTree.end())
                {
                    lineSegmentIt = lineSegmentTree.insert(lineSegment).first;
                }
                else if (cp > 0 && lineSegmentIt != lineSegmentTree.end())
                {
                    lineSegmentTree.erase(lineSegmentIt);
                    lineSegmentIt = lineSegmentTree.end();
                }
            };

            processLineSegment(cwLineSegment);
            processLineSegment(ccwLineSegment);
        }

        return visiblePoints;
    }

    static bool areVisible(const Point &p1, const Point &p2, const vec<Polygon> &obstacles) {
        LineSegment visionLine(p1, p2);
        bool visible = true;

        for (auto it = obstacles.begin(); visible && it != obstacles.end(); ++it)
        {
            const Polygon &polygon = *it;
            for (auto pointIt = polygon.begin(); pointIt != polygon.end(); ++pointIt)
            {
                auto otherPointIt = pointIt + 1 == polygon.end() ? polygon.begin() : pointIt + 1;
                LineSegment edge(*pointIt, *otherPointIt);

                Point p, q;

                if (visionLine.intersection(edge, p, q) == LineSegment::IntersectionMode::PointIntersection)
                {
                    bool pointIntersection = approx(p.x(), pointIt->x()) && approx(p.y(), pointIt->y());
                    bool otherPointIntersection = approx(p.x(), otherPointIt->x()) && approx(p.y(), otherPointIt->y());

                    visible = pointIntersection || otherPointIntersection;
                }
            }
        }

        return visible;
    }

    static void addEdgesToGraph(Graph<Point> &graph, const Point &vertex, const vec<Point> &vertices, bool mirror) {
        for (const auto & vertice : vertices) {
            graph.addEdge(vertex, vertice);
            if (mirror) graph.addEdge(vertice, vertex);
        }
    }

    struct PointRef
    {
        size_t polygonId;
        size_t localPointId;
        size_t globalPointId;
    };

    struct LineSegmentRef
    {
        size_t lineSegmentId;
        size_t polarPointAId;
        size_t polarPointBId;
    };

    static Graph<Point> find(const Point &start, const Point &destination, const vec<Polygon> &obstacles) {
        Graph<Point> graph;

        PointRef externPoint = { obstacles.size(), 0, 0 };

        addEdgesToGraph(graph, start, visibleVertices(start, obstacles, externPoint), true);
        addEdgesToGraph(graph, destination, visibleVertices(destination, obstacles, externPoint), true);

        size_t pointCounter = 0;

        for (auto it = obstacles.begin(); it != obstacles.end(); ++it)
        {
            const Polygon &polygon = *it;
            for (auto point = polygon.begin(); point != polygon.end(); ++point)
            {
                PointRef ref =
                        {
                                static_cast<size_t>(std::distance(obstacles.begin(), it)),
                                static_cast<size_t>(std::distance(polygon.begin(), point)),
                                pointCounter++
                        };
                addEdgesToGraph(graph, *point, visibleVertices(*point, obstacles, ref), true);
            }
        }

        if (areVisible(start, destination, obstacles))
        {
            graph.addEdge(start, destination);
            graph.addEdge(destination, start);
        }

        return graph;
    }
};

struct ShortestPath {
    struct SPNode
    {
        ftype distance;
        Point point;
        Point parent;
    };

    static Path find(const Point& start, const Point& destination, const Graph<Point> &visibilityGraph) {
        auto comp = [](const SPNode &a, const SPNode &b)
        {
            if (a.distance != b.distance)
                return a.distance > b.distance;

            return a.point < b.point;
        };

        typedef std::priority_queue<SPNode, std::vector<SPNode>, decltype(comp)> PointQueue;
        PointQueue queue(comp, std::vector<SPNode>());

        queue.push(SPNode {0.0, start, Point(0.0, 0.0)});

        std::map<Point, SPNode> sp;

        while (!queue.empty())
        {
            SPNode top = queue.top();
            queue.pop();

            if (sp.count(top.point))
                continue;

            sp[top.point] = top;

            const auto &edges = visibilityGraph.getEdges(top.point);

            for (const auto & edge : edges) {
                Point diff = top.point - edge;

                if (!sp.count(edge))
                    queue.push({ top.distance + std::hypot(diff.x(), diff.y()), edge, top.point });
            }
        }

        Path path;
        Point d = destination;

        while (d != start)
        {
            path.push_back(d);
            d = sp.at(d).parent;
        }

        path.push_back(start);
        std::reverse(path.begin(), path.end());

        return path;
    }
};


void solve() {
    int n;
    cin >> n;
    Polygon a(n);
    for (auto& p: a) cin >> p;
    Point g, t;
    cin >> g >> t;
    const auto& vis = VisibilityGraph::find(g, t, {a});
    for (auto& p: vis) {
        cout << p.F << "------\n";
        for (auto& v: p.S) {
            cout << v << nl;
        }
    }
//    const auto& path = ShortestPath::find(g, t, vis);
//    for (auto& p: path) {
//        cout << "(" << p.x() << ", " << p.y() << ")" << nl;
//    }
//    cout << "DONE" << nl;
}

int32_t main() {
    fast
    freopen("./in", "r", stdin);
    solve();
}
