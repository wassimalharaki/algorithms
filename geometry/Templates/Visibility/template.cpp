/**
 * 16:01:05 8/26/24
 * template
 */
// ./ICPC/Geometry/Templates/template.cpp
#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <utility>

using namespace std;
using namespace __gnu_pbds;
#define int long long
#define uint unsigned int
#define double long double
#define fast ios_base::sync_with_stdio(false);cin.tie(NULL);
#define nl '\n'
#define NO (cout << "NO" << nl);
#define YES (cout << "YES" << nl);
#define F first
#define S second
#define INF 1000000000000000000ll
#define MOD 1000000007ll
#define pii pair<int, int>
#define P complex<int>
#define X real()
#define Y imag()
#define vec vector
#define LT int T; cin >> T; while (T--)
#define ftype double
#define ai2 array<int, 2>
template<typename T>
using vec2d = vector<vector<T>>;
template<typename T>
using ordered_set = tree<T, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;
using indexed_set = tree<int, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;

// if the sign of the area is negative, the points turn clockwise.
//P centroid(const vec<P>& a) {
//    P r{0, 0};
//    for (auto& p: a) r += p;
//    return r / a.size();
//}
// Alternative equation for center of gravity
// Cx = 1/6A [sum from i = 0; i < n - 1]((x[i] + x[i + 1]) * (x[i] * y[i + 1] - x[i + 1] * y[i]))
// Cy = 1/6A [sum from i = 0; i < n - 1]((y[i] + y[i + 1]) * (x[i] * y[i + 1] - x[i + 1] * y[i]))
// Where A is the area of the polygon.

enum class Orientation {
    left_turn = 1,
    right_turn = -1,
    collinear = 0
};
const double eps = 1e-12l;

bool equal(double a, double b) {
    return abs(a - b) <= max(abs(a), abs(b)) * eps;
}

bool lessThan(double a, double b) {
    return (b - a) > max(abs(a), abs(b)) * eps;
}

struct point2d {
    ftype x, y;
    point2d(): x{0}, y{0} {}
    point2d(ftype x, ftype y): x(x), y(y) {}
    bool operator==(const point2d& p) const {
        return equal(x, p.x) and equal(y, p.y);
    }
    bool operator<(const point2d& p) const {
        return equal(x, p.x) ? lessThan(y, p.y) : lessThan(x, p.x);
    }
    point2d& operator+=(const point2d &t) {
        x += t.x;
        y += t.y;
        return *this;
    }
    point2d& operator-=(const point2d &t) {
        x -= t.x;
        y -= t.y;
        return *this;
    }
    point2d& operator*=(ftype t) {
        x *= t;
        y *= t;
        return *this;
    }
    point2d& operator/=(ftype t) {
        x /= t;
        y /= t;
        return *this;
    }
    point2d operator+(const point2d &t) const {
        return point2d(*this) += t;
    }
    point2d operator-(const point2d &t) const {
        return point2d(*this) -= t;
    }
    point2d operator*(ftype t) const {
        return point2d(*this) *= t;
    }
    double operator*(const point2d& p) const { // Cross product
        return x * p.y - y * p.x;
    }
    point2d operator/(ftype t) const {
        return point2d(*this) /= t;
    }
    [[nodiscard]] double dot(const point2d& p) const {
        return x * p.x + y * p.y;
    }
    [[nodiscard]] double norm() const {
        return dot(*this);
    }
    [[nodiscard]] double abs() const {
        return sqrt(norm());
    }
    [[nodiscard]] point2d normalize() const {
        double v = abs();
        if (::abs(v) < eps) return *this;
        return *this / v;
    }
    [[nodiscard]] point2d normal() const {
        return { -y, x };
    }
};
point2d operator*(ftype a, point2d b) {
    return b * a;
}
istream& operator>>(istream& in, point2d& p) {
    int x, y;
    in >> x >> y;
    p.x = x;
    p.y = y;
    return in;
}
ostream& operator<<(ostream& out, const point2d& p) {
    return out << "(" << p.x << ", " << p.y << ")";
}
struct point3d {
    ftype x, y, z;
    point3d() {}
    point3d(ftype x, ftype y, ftype z): x(x), y(y), z(z) {}
    point3d& operator+=(const point3d &t) {
        x += t.x;
        y += t.y;
        z += t.z;
        return *this;
    }
    point3d& operator-=(const point3d &t) {
        x -= t.x;
        y -= t.y;
        z -= t.z;
        return *this;
    }
    point3d& operator*=(ftype t) {
        x *= t;
        y *= t;
        z *= t;
        return *this;
    }
    point3d& operator/=(ftype t) {
        x /= t;
        y /= t;
        z /= t;
        return *this;
    }
    point3d operator+(const point3d &t) const {
        return point3d(*this) += t;
    }
    point3d operator-(const point3d &t) const {
        return point3d(*this) -= t;
    }
    point3d operator*(ftype t) const {
        return point3d(*this) *= t;
    }
    point3d operator/(ftype t) const {
        return point3d(*this) /= t;
    }
};
point3d operator*(ftype a, point3d b) {
    return b * a;
}
using PolySegment = pair<point2d, point2d>;
struct Poly {
    vec<PolySegment> s;
};
struct PolyP {
    vec<point2d> p;
};
//double ternary_search(double l, double r) {
//    double eps = 1e-9;              //set the error limit here
//    while (r - l > eps) {
//        double m1 = l + (r - l) / 3;
//        double m2 = r - (r - l) / 3;
//        double f1 = f(m1);      //evaluates the function at m1
//        double f2 = f(m2);      //evaluates the function at m2
//        if (f1 < f2)
//            l = m1;
//        else
//            r = m2;
//    }
//    return f(l);                    //return the maximum of f(x) in [l, r]
//}


Orientation getOrientation(const point2d& a, const point2d& b, const point2d& c) {
    auto det = (b - a) * (c - a);
    return static_cast<Orientation>(
            lessThan(0, det) - lessThan(det, 0)
    );
}

struct Ray {
    point2d org, dir;
    Ray(): org{}, dir{} {};
    Ray(point2d orig, point2d dir): org{orig}, dir{dir} {}

    /** Find the nearest intersection point of ray and line segment.
         * @param segment
         * @param out_point reference to a variable where the nearest
         *        intersection point will be stored (can be changed even
         *        when there is no intersection)
         * @return true iff the ray intersects the line segment
         */
    bool intersects(const PolySegment& seg, point2d& outp) const {
        auto ao = org - seg.F;
        auto ab = seg.S - seg.F;
        auto det = ab * dir;
        if (equal(det, 0)) {
            auto abo = getOrientation(seg.F, seg.S, org);
            if (abo != Orientation::collinear) return false;
            auto dist_a = ao.dot(dir);
            auto dist_b = (org - seg.S).dot(dir);

            if (dist_a > 0 and dist_b > 0) return false;
            else if ((dist_a > 0) != (dist_b > 0)) outp = org;
            else if (dist_a > dist_b) outp = seg.F;
            else outp = seg.S;
            return true;
        }
        auto u  = (ao * dir) / det;
        if (lessThan(u, 0) or lessThan(1, u)) return false;
        auto t = -(ab * ao) / det;
        outp = org + t * dir;
        return equal(t, 0) or t > 0;
    }
};

ftype dot(point2d a, point2d b) {
    return a.dot(b);
}
ftype dot(point3d a, point3d b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
// Norm(a) = abs(a * a) = dot(a, a)

// Projection of a onto b: dot(a, b) / len(b)

// Angle between vectors = arccos(dot(a, b) / (len(a) * len(b)))

// From the previous point we may see that the dot product is positive if the angle between them is acute,
// negative if it is obtuse and it equals zero if they are orthogonal, i.e. they form a right angle.
ftype norm(point2d a) {
    return a.norm();
}
double abs(point2d a) {
    return a.abs();
}
double proj(point2d a, point2d b) {
    return dot(a, b) / abs(b);
}
//double angle(const P& p1, const P& p0, const P& p2) {
//    double p01 = dist(p0, p1), p02 = dist(p0, p2), p12 = dist(p1, p2);
//    return acos((sqr(p01) + sqr(p02) - sqr(p12)) / (2 * p01 * p02));
//}
//
//double minorAngle(const seg& s, const P& p) {
//    double r0 = angle({0, 0}, p, s.F < s.S ? s.F : s.S);
//    return r0;
//}

double angle(point2d a, point2d b) {
    return atan2(b.y - a.y, b.x - a.x);
}
// To see the next important property we should take a look at the set of points r
// for which dot(r, a) = C for some fixed constant C.
// You can see that this set of points is exactly the set of points for which the projection onto
// a is the point dot(C, a / len(a)) and they form a hyperplane orthogonal to a.
// In 2D these vectors will form a line, in 3D they will form a plane.
// Note that this result allows us to define a line in 2D as dot(r, n) = C
// or dot(r - r0, n) = 0. Where n is vector orthogonal to the line and r0 is any vector already present on the line and C = dot(r0, n).
// In the same manner a plane can be defined in 3D.
point3d cross(point3d a, point3d b) {
    return point3d(a.y * b.z - a.z * b.y,
                   a.z * b.x - a.x * b.z,
                   a.x * b.y - a.y * b.x);
}
ftype triple(point3d a, point3d b, point3d c) {
    return dot(a, cross(b, c));
}
ftype cross(point2d a, point2d b) {
    return a * b;
}
point2d intersect(point2d a1, point2d d1, point2d a2, point2d d2) {
    return a1 + cross(a2 - a1, d2) / cross(d1, d2) * d1;
}
point3d intersect(point3d a1, point3d n1, point3d a2, point3d n2, point3d a3, point3d n3) {
    point3d x(n1.x, n2.x, n3.x);
    point3d y(n1.y, n2.y, n3.y);
    point3d z(n1.z, n2.z, n3.z);
    point3d d(dot(a1, n1), dot(a2, n2), dot(a3, n3));
    return point3d(triple(d, y, z),
                   triple(x, d, z),
                   triple(x, y, d)) / triple(n1, n2, n3);
}
int signed_area_parallelogram(point2d p1, point2d p2, point2d p3) {
    return cross(p2 - p1, p3 - p2);
}

double triangle_area(point2d p1, point2d p2, point2d p3) {
    return abs(signed_area_parallelogram(p1, p2, p3)) / 2.0;
}

bool clockwise(point2d p1, point2d p2, point2d p3) {
    return signed_area_parallelogram(p1, p2, p3) < 0;
}

bool counter_clockwise(point2d p1, point2d p2, point2d p3) {
    return signed_area_parallelogram(p1, p2, p3) > 0;
}

int pythagoras(int a = 0, int b = 0, int c = 0) {
    if (a == 0) {
        return c * c - b * b;
    } else if (b == 0) {
        return c * c - a * a;
    } else {
        return a * a + b * b;
    }
}

double pythagoras(double a = 0, double b = 0, double c = 0) {
    if (a == 0) {
        return sqrt(c * c - b * b);
    } else if (b == 0) {
        return sqrt(c * c - a * a);
    } else {
        return sqrt(a * a + b * b);
    }
}

template<class T>
T sqr(T v) { return v * v; }

double PPDistance(point2d p0, point2d p1) {
    return sqrt(sqr(p0.x - p1.x) + sqr(p0.y - p1.y));
}

double PLDistance(point2d l0, point2d l1, point2d p) {
    double x0 = p.x, y0 = p.y,
           x1 = l0.x, y1 = l0.y,
           x2 = l1.x, y2 = l1.y;
    return abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / sqrt(sqr(y2 - y1) + sqr(x2 - x1));
}

double PSDistance(point2d s0, point2d s1, point2d p) {
    double x = p.x, y = p.y,
           x1 = s0.x, y1 = s0.y,
           x2 = s1.x, y2 = s1.y;
    double A = x - x1;
    double B = y - y1;
    double C = x2 - x1;
    double D = y2 - y1;

    double dot = A * C + B * D;
    double len_sq = C * C + D * D;
    double param = -1;
    if (len_sq != 0) param = dot / len_sq;

    double xx, yy;
    if (param < 0) {
        xx = x1;
        yy = y1;
    } else if (param > 1) {
        xx = x2;
        yy = y2;
    } else {
        xx = x1 + param * C;
        yy = y1 + param * D;
    }

    double dx = x - xx;
    double dy = y - yy;
    return sqrt(dx * dx + dy * dy);
}
PolySegment PSDistanceVec(point2d s0, point2d s1, point2d p) {
    double x = p.x, y = p.y,
            x1 = s0.x, y1 = s0.y,
            x2 = s1.x, y2 = s1.y;
    double A = x - x1;
    double B = y - y1;
    double C = x2 - x1;
    double D = y2 - y1;

    double dot = A * C + B * D;
    double len_sq = C * C + D * D;
    double param = -1;
    if (len_sq != 0) param = dot / len_sq;

    double xx, yy;
    if (param < 0) {
        xx = x1;
        yy = y1;
    } else if (param > 1) {
        xx = x2;
        yy = y2;
    } else {
        xx = x1 + param * C;
        yy = y1 + param * D;
    }
    return {p, {xx, yy}};
}

double PLDistance(point2d p, double a, double b, double c) {
    double x0 = p.x, y0 = p.y;
    return abs(a * x0 + b * y0 + c) / sqrt(sqr(a) + sqr(b));
}

point2d ClosestPointToLine(point2d p, double a, double b, double c) {
    double x0 = p.x, y0 = p.y;
    double x = (b * (b * x0 - a * y0) - a * c) / (sqr(a) + sqr(b));
    double y = (a * (-b * x0 + a * y0) - b * c) / (sqr(a) + sqr(b));
    return {x, y};
}

// In case they are parallel and the lines have the form ax + by + c = 0
double LLDistance(double c1, double c2, double a1, double b1) {
    return abs(c1 - c2) / sqrt(sqr(a1) + sqr(b1));
}

// Returns the equation of the line in the form of y = ax + b
// First element is a, second is b
pair<double, double> Vec2Norm(point2d p0, point2d p1) {
    double a = (p0.y - p1.y) / (p0.x - p1.x);
    // using: y - y1 = a(x - x1), calculate b
    // y = ax - ax1 + y1; b = y1 - a.x1
    double b = p0.y - a * p0.x;
    return {a, b};
}

// Given three collinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(point2d p, point2d q, point2d r)
{
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
        return true;

    return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(point2d p, point2d q, point2d r)
{
    // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
    // for details of below formula.
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);

    if (val == 0) return 0;  // collinear

    return (val > 0)? 1: 2; // clock or counterclock wise
}

double slope(point2d p1, point2d q1) {
    return (p1.y - q1.y) / (p1.x - q1.x);
}

// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool doIntersect(point2d p1, point2d q1, point2d p2, point2d q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4) {
        return true;
    }

    // Special Cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false; // Doesn't fall in any of the above cases
}

bool doIntersect(const PolySegment& s, const Poly& p) {
    for (auto& seg: p.s) {
        if (doIntersect(seg.F, seg.S, s.F, s.S)) return true;
    }
    return false;
}

// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines
// intersect the intersection point may be stored in second element of the pair.
pair<char, point2d> segmentSlopeIntersection(point2d p0, point2d p1, point2d p2, point2d p3) {
    double p0_x = p0.x, p0_y = p0.y;
    double p1_x = p1.x, p1_y = p1.y;
    double p2_x = p2.x, p2_y = p2.y;
    double p3_x = p3.x, p3_y = p3.y;

    double s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    double s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
        // Collision detected
        return {1, {p0_x + (t * s1_x), p0_y + (t * s1_y)}};
    }

    return {0, {}}; // No collision
}
struct Intersection {
    bool fail = false;
    double s0{}, s1{};
    point2d p{};
};

Intersection segmentProperIntersection(const pair<point2d, point2d>& s0, const pair<point2d, point2d>& s1) {
    point2d p0 = s0.F, p1 = s0.S, p2 = s1.F, p3 = s1.S;
    double p0_x = p0.x, p0_y = p0.y;
    double p1_x = p1.x, p1_y = p1.y;
    double p2_x = p2.x, p2_y = p2.y;
    double p3_x = p3.x, p3_y = p3.y;

    double s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    double s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (isinf(s) or isinf(t)) {
        return {true};
    }

    return {false, t, s, {p0_x + (t * s1_x), p0_y + (t * s1_y)}};
}


double circleArea(double r, double a = 2 * M_PI) {
    return M_PI * r * r * (a / (2 * M_PI));
}
// Heron's formula
double triangle(double a, double b, double c) {
    double s = (a + b + c) / 2;
    return sqrt(s * (s - a) * (s - b) * (s - c));
}
double squareArea(double a) { return sqr(a); }

double segmentDistance(point2d p1, point2d q1,
                       point2d p2, point2d q2) {
    vec<double> tmp {
            PSDistance(p1, q1, p2),
            PSDistance(p1, q1, q2),
            PSDistance(p2, q2, p1),
            PSDistance(p2, q2, q1),
    };
    return *min_element(tmp.begin(), tmp.end());
}

double segmentDistance(const PolySegment& s0, const PolySegment& s1) {
    return segmentDistance(s0.F, s0.S, s1.F, s1.S);
}


double PolyPDistance(point2d x, const Poly& poly) {
    double mn = INF;
    for (auto& seg: poly.s) {
        mn = min(mn, PSDistance(seg.F, seg.S, x));
//        auto &p1 = seg.p0,
//             &p2 = seg.p1;
//        double r = dot(p1, p2);
//        r /= norm(p2 - p1);
//        double dist;
//        if (r < 0) {
//            dist = abs(x - p1);
//        } else if (r > 1) {
//            dist = abs(p2 - x);
//        } else {
//            dist = sqrt(norm(x - p1) - sqr(r * abs(p2 - p1)));
//        }
//        mn = min(mn, dist);
    }
    return mn;
}
double len(const PolySegment& s) {
    return PPDistance(s.F, s.S);
}
PolySegment PolyPDistanceVec(point2d x, const Poly& poly) {
    PolySegment mn{{0, 0}, {INF, INF}};
    for (auto& seg: poly.s) {
        const auto& v = PSDistanceVec(seg.F, seg.S, x);
        if (len(v) < len(mn)) mn = v;
    }
    return mn;
}

bool isin(point2d x, const PolyP& poly) {
    double minX = INF, minY = INF, maxX = -INF, maxY = -INF;
    for (auto& p: poly.p) {
        minX = min(minX, p.x);
        minY = min(minY, p.y);
        maxX = max(maxX, p.x);
        maxY = max(maxY, p.y);
    }
    if (x.x < minX or x.x > maxX or x.y < minY or x.y > maxY) return false;
    bool in = false;
    int sz = (int)poly.p.size();
    for (int i = 0, j = sz - 1; i < sz; j = i++) {
        if (
                ((poly.p[i].y > x.y) ^ (poly.p[j].y > x.y))
                and (x.x < ( poly.p[j].x - poly.p[i].x ) * ( x.y - poly.p[i].y ) / ( poly.p[j].y - poly.p[i].y ) + poly.p[i].x)
                ) {
            in = not in;
        }
    }
    return in;
}
bool ison(point2d x, const Poly& poly) {
    return ranges::any_of(
            poly.s, [&](auto& s) {
        return onSegment(s.F, x, s.S);
    });
}

double PolyPolyDistance(const Poly& p0, const Poly& p1) {
    double mn = INF;
    for (auto& s0: p0.s) {
        for (auto& s1: p1.s) {
            mn = min(mn, segmentDistance(s0, s1));
        }
    }
    return mn;
}
point2d rotate(const point2d& p0, double a, const point2d& ref = {0, 0}) {
    point2d res;
    double xs = p0.x - ref.x, ys = p0.y - ref.y;
    res.x = ref.x + (xs * cos(a) - ys * sin(a));
    res.y = ref.y + (xs * sin(a) + ys * cos(a));
    return res;
}

PolySegment rotate(const PolySegment& s0, double a, const point2d& ref = {0, 0}) {
    return {rotate(s0.F, a, ref), rotate(s0.S, a, ref)};
}

Poly rotate(const Poly& p, double a, const point2d& ref = {0, 0}) {
    Poly res;
    for (auto& s: p.s) {
        res.s.push_back(rotate(s, a, ref));
    }
    return res;
}

Poly convert(const PolyP& p) {
    Poly res;
    int n = (int)p.p.size();
    for (int i = 0; i < n; i++) {
        res.s.emplace_back(p.p[i], p.p[(i + 1) % n]);
    }
    return res;
}
PolyP convert(const Poly& p) {
    PolyP res;
    int n = (int)p.s.size();
    for (int i = 0; i < n; i++) {
        res.p.emplace_back(p.s[i].F);
    }
    return res;
}

/*
 * TODO
 * https://byjus.com/maths/area-segment-circle/
 * https://www.geometrictools.com/Documentation/DistanceLine3Line3.pdf
 * https://en.wikipedia.org/wiki/Heronian_triangle
 * area of a section of an ellipse formula
 * https://byjus.com/maths/area-of-ellipse/
 */
// to go from general form: x^2 + y^2 + Ax + By + C = 0
// to standard form: (x−a)^2 + (y−b)^2 = r^2
// Complete the sequare
/*
 * Equation of a circle (x - x_r) ** 2 + (y - y_r) ** 2 = r ** 2
 * A(2pi) = pi * r ^ 2
 * C = 2 * pi * r
 * A(a) = A(2pi) * (a / 2pi)
 *
 * Start with:x2 + y2 − 2x − 4y − 4 = 0
 * Put xs and ys together:(x2 − 2x) + (y2 − 4y) − 4 = 0
 * Constant on right:(x2 − 2x) + (y2 − 4y) = 4
 *
 * Now complete the square for x (take half of the −2, square it, and add to both sides):
 *
 * (x2 − 2x + (−1)2) + (y2 − 4y) = 4 + (−1)2
 * And complete the square for y (take half of the −4, square it, and add to both sides):
 * (x2 − 2x + (−1)2) + (y2 − 4y + (−2)2) = 4 + (−1)2 + (−2)2
 */


using seg = pair<point2d, point2d>;
int sign(double v) {
    if (equal(v, 0)) return 0;
    else if (lessThan(v, 0)) return -1;
    else return 1;
}

enum Side {
    Right = -1,
    On = 0,
    Left = 1
};
Side side(const seg& s, const point2d& p) {
    double Bx = s.S.x, Ax = s.F.x, By = s.S.y, Ay = s.F.y;
    return (Side)sign((Bx - Ax) * (p.y - Ay) - (By - Ay) * (p.x - Ax));
}

point2d segP(const seg& s, double t) {
    double x1 = s.F.x, y1 = s.F.y, x2 = s.S.x, y2 = s.S.y;
    return {x1 + t * (x2 - x1), y1 + t * (y2 - y1)};
}

// Equation of a vector (x1, y1), (x2, y2)
// is: x(t) = x1 + t(x2 - x1), y(t) = y1 + t(y2 - y1)
vec<point2d> intersection(PolySegment& s, double r, point2d c) {
    double x1 = s.F.x, y1 = s.F.y, x2 = s.S.x, y2 = s.S.y, xc = c.x, yc = c.y;
    auto gx = [&](double t) {
        return x1 + t * (x2 - x1);
    };
    auto gy = [&](double t) {
        return y1 + t * (y2 - y1);
    };
    double A = sqr(x2 - x1) + sqr(y2 - y1),
           B = 2 * ((x2 - x1) * (x1 - xc) + (y2 - y1) * (y1 - yc)),
           C = sqr(x1 - xc) + sqr(y1 - yc) - r * r;
    double delta = B * B - 4 * A * C;
    if (delta < 0) return {};
    else if (delta == 0) {
        double t = -B / (2 * A);
        return {{ gx(t), gy(t) }};
    }
    delta = sqrt(delta);
    vec<point2d> res;
    double t0 = (-B + delta) / (2 * A), t1 = (-B - delta) / (2 * A);
    if (0 <= t0 and t0 <= 1) res.emplace_back(gx(t0), gy(t0));
    if (0 <= t1 and t1 <= 1) res.emplace_back(gx(t1), gy(t1));
    return res;
}

double angle(point2d p0, point2d p1, point2d p2) {
    double a1 = angle(p0, p1), a2 = angle(p1, p2);
    double a3 = a1 - a2;
    if (a3 < 0) a3 += 2 * M_PI;
    if (a3 > 2 * M_PI) a3 -= 2 * M_PI;
    return a3;
}
double toDeg(double a) { return a / M_PI * 180; }


struct LineSegmentDistComparer {
    using segment_type = PolySegment;
    point2d origin;

    explicit LineSegmentDistComparer(point2d origin) :
    origin(origin)
            {
            }

    /** Check whether the line segment x is closer to the origin than the
     * line segment y.
     * @param x line segment: left hand side of the comparison operator
     * @param y line segment: right hand side of the comparison operator
     * @return true iff x < y (x is closer than y)
     */
    bool operator()(const segment_type& x, const segment_type& y) const {
        auto a = x.F, b = x.S;
        auto c = y.F, d = y.S;

        assert(
                getOrientation(origin, a, b) != Orientation::collinear &&
                "AB must not be collinear with the origin.");
        assert(
                getOrientation(origin, c, d) != Orientation::collinear &&
                "CD must not be collinear with the origin.");

        // sort the endpoints so that if there are common endpoints,
        // it will be a and c
        if (b == c or b == d)
            std::swap(a, b);
        if (a == d)
            std::swap(c, d);

        // cases with common endpoints
        if (a == c)
        {
            auto oad = getOrientation(origin, a, d);
            auto oab = getOrientation(origin, a, b);
            if (b == d || oad != oab)
                return false;
            return getOrientation(a, b, d) != getOrientation(a, b, origin);
        }

        // cases without common endpoints
        auto cda = getOrientation(c, d, a);
        auto cdb = getOrientation(c, d, b);
        if (cdb == Orientation::collinear && cda == Orientation::collinear) {
            return (origin - a).norm() < (origin - c).norm();
        }
        else if (cda == cdb ||
                 cda == Orientation::collinear ||
                 cdb == Orientation::collinear)
        {
            auto cdo = getOrientation(c, d, origin);
            return cdo == cda || cdo == cdb;
        }
        else
        {
            auto abo = getOrientation(a, b, origin);
            return abo != getOrientation(a, b, c);
        }
    }
};
struct AngleComparer {
    point2d vertex;

    explicit AngleComparer(point2d origin) : vertex(origin) {}

    bool operator()(const point2d& a, const point2d& b) const
    {
        auto is_a_left = lessThan(a.x, vertex.x);
        auto is_b_left = lessThan(b.x, vertex.x);
        if (is_a_left != is_b_left)
            return is_b_left;

        if (equal(a.x, vertex.x) && equal(b.x, vertex.x))
        {
            if (!lessThan(a.y, vertex.y) ||
                !lessThan(b.y, vertex.y))
            {
                return lessThan(b.y, a.y);
            }
            return lessThan(a.y, b.y);
        }

        auto oa = a - vertex;
        auto ob = b - vertex;
        auto det = cross(oa, ob);
        if (equal(det, 0.f))
        {
            return oa.norm() < ob.norm();
        }
        return det < 0;
    }
};
struct VisibilityEvent {
    enum EventType {
        start_vertex, end_vertex
    };
    EventType type;
    PolySegment segment;
    VisibilityEvent() = default;
    VisibilityEvent(EventType type, PolySegment segment): type(type), segment(std::move(segment)) {}
    [[nodiscard]] const auto& point() const { return segment.F; }
};

PolyP visibility_polygon(point2d point, const Poly& p) {
    using segment_type = PolySegment;
    using event_type = VisibilityEvent;
    using segment_comparer_type = LineSegmentDistComparer;
    segment_comparer_type cmp_dist{ point };
    std::set<segment_type, segment_comparer_type> state{ cmp_dist };
    std::vector<event_type> events;
    for (auto& segment: p.s)
    {
        // Sort line segment endpoints and add them as events
        // Skip line segments collinear with the point
        auto pab = getOrientation(point, segment.F, segment.S);
        if (pab == Orientation::collinear)
        {
            continue;
        }
        else if (pab == Orientation::right_turn)
        {
            events.emplace_back(event_type::start_vertex, segment);
            events.emplace_back(
                    event_type::end_vertex,
                    segment_type{ segment.S, segment.F });
        }
        else
        {
            events.emplace_back(
                    event_type::start_vertex,
                    segment_type{ segment.S, segment.F });
            events.emplace_back(event_type::end_vertex, segment);
        }

        // Initialize state by adding line segments that are intersected
        // by vertical ray from the point
        auto a = segment.F, b = segment.S;
        if (a.x > b.x)
            std::swap(a, b);

        auto abp = getOrientation(a, b, point);
        if (abp == Orientation::right_turn &&
            (equal(b.x, point.x) ||
             (a.x < point.x && point.x < b.x)))
        {
            state.insert(segment);
        }
    }

    // sort events by angle
    AngleComparer cmp_angle{ point };
    std::sort(events.begin(), events.end(), [&cmp_angle](auto&& a, auto&& b)
    {
        // if the points are equal, sort end vertices first
        if (a.point() == b.point())
            return a.type == event_type::end_vertex &&
                   b.type == event_type::start_vertex;
        return cmp_angle(a.point(), b.point());
    });

    // find the visibility polygon
    vec<point2d> vertices;
    for (auto&& event : events)
    {
        if (event.type == event_type::end_vertex)
            state.erase(event.segment);

        if (state.empty())
        {
            vertices.push_back(event.point());
        }
        else if (cmp_dist(event.segment, *state.begin()))
        {
            // Nearest line segment has changed
            // Compute the intersection point with this segment
            point2d intersection;
            Ray ray{ point, event.point() - point };
            auto nearest_segment = *state.begin();
            auto intersects = ray.intersects(nearest_segment, intersection);
            assert(intersects &&
                           "Ray intersects line segment L iff L is in the state");

            if (event.type == event_type::start_vertex)
            {
                vertices.push_back(intersection);
                vertices.push_back(event.point());
            }
            else
            {
                vertices.push_back(event.point());
                vertices.push_back(intersection);
            }
        }

        if (event.type == event_type::start_vertex)
            state.insert(event.segment);
    }

    // remove collinear points
    auto top = vertices.begin();
    for (auto it = vertices.begin(); it != vertices.end(); ++it)
    {
        auto prev = top == vertices.begin() ? vertices.end() - 1 : top - 1;
        auto next = it + 1 == vertices.end() ? vertices.begin() : it + 1;
        if (getOrientation(*prev, *it, *next) != Orientation::collinear)
            *top++ = *it;
    }
    vertices.erase(top, vertices.end());
    for (auto& s: p.s) {
        if (s.F == point or s.S == point) {
            vertices.push_back(point);
            break;
        }
    }
    return {vertices};
}

double mnAns = 1e20l;
pair<point2d, point2d> mnAnsP;
void updateMnAns(const point2d& p0, const point2d& p1) {
    double r = (p1 - p0).abs();
    if (lessThan(r, mnAns)) {
        mnAnsP = {p0, p1};
        mnAns = r;
    }
}
vec<point2d> t;

// Returns the closes pair of points
void minDist(vec<point2d>& a, int l = 0, int r = -1) {
    if (r == -1) {
        r = (int)a.size();
        t.resize(a.size());
        sort(a.begin(), a.end());
    }
    if (r - l <= 3) {
        for (int i = l; i < r; i++) {
            for (int j = i + 1; j < r; j++) {
                updateMnAns(a[i], a[j]);
            }
        }
        sort(a.begin() + l, a.begin() + r, [&](auto& p0, auto& p1) {
            return lessThan(p0.y, p1.y);
        });
        return;
    }
    int m = (l + r) >> 1;
    double midx = a[m].x;
    minDist(a, l, m);
    minDist(a, m, r);
    merge(a.begin() + l, a.begin() + m, a.begin() + m, a.begin() + r, t.begin(), [&](auto& p0, auto& p1) {
        return lessThan(p0.y, p1.y);
    });
    copy(t.begin(), t.begin() + r - l, a.begin() + l);

    int tsz = 0;
    for (int i = l; i < r; i++) {
        if (lessThan(a[i].x - midx, mnAns)) {
            for (int j = tsz - 1; j >= 0 && lessThan(a[i].y - t[j].y, mnAns); --j)
                updateMnAns(a[i], t[j]);
            t[tsz++] = a[i];
        }
    }
}

// TEST CASES
// _____________________________________________________
namespace TestCases {
    void REQUIRE(bool v) {
        if (not v) {
            cerr << "INVALID EXPECTED TRUE GOT FALSE" << nl;
            exit(1);
        }
    }

    void REQUIRE_FALSE(bool v) {
        if (v) {
            cerr << "INVALID EXPECTED FALSE GOT TRUE" << nl;
            exit(1);
        }
    }

    bool equal(const point2d &a, const point2d &b) {
        return a == b;
    }

    using vector_type = point2d;
    using segment_type = PolySegment;
    using segment_comparer_type = LineSegmentDistComparer;
    using angle_comparer_type = AngleComparer;

    void test_line_segment_is_closer(
            const segment_comparer_type &cmp,
            vector_type a, vector_type b,
            vector_type c, vector_type d) {
        REQUIRE(cmp(segment_type{a, b}, {c, d}));
        REQUIRE(cmp(segment_type{b, a}, {c, d}));
        REQUIRE(cmp(segment_type{a, b}, {d, c}));
        REQUIRE(cmp(segment_type{b, a}, {d, c}));

        REQUIRE_FALSE(cmp(segment_type{c, d}, {a, b}));
        REQUIRE_FALSE(cmp(segment_type{d, c}, {a, b}));
        REQUIRE_FALSE(cmp(segment_type{c, d}, {b, a}));
        REQUIRE_FALSE(cmp(segment_type{d, c}, {b, a}));
    }

    void test_line_segments_are_equal(
            const segment_comparer_type &cmp,
            vector_type a, vector_type b,
            vector_type c, vector_type d) {
        REQUIRE_FALSE(cmp(segment_type{a, b}, {c, d}));
        REQUIRE_FALSE(cmp(segment_type{b, a}, {c, d}));
        REQUIRE_FALSE(cmp(segment_type{a, b}, {d, c}));
        REQUIRE_FALSE(cmp(segment_type{b, a}, {d, c}));

        REQUIRE_FALSE(cmp(segment_type{c, d}, {a, b}));
        REQUIRE_FALSE(cmp(segment_type{d, c}, {a, b}));
        REQUIRE_FALSE(cmp(segment_type{c, d}, {b, a}));
        REQUIRE_FALSE(cmp(segment_type{d, c}, {b, a}));
    }

    void T0() {
        cout << "Compare 2 line segments with no common endpoints" << nl;
        segment_comparer_type cmp{{0, 0}};

        test_line_segment_is_closer(cmp, {1, 1}, {1, -1}, {2, 1}, {2, -1});
        test_line_segment_is_closer(cmp, {1, 1}, {1, -1}, {2, 2}, {2, 3});
    }

    void T1() {
        cout << "Compare 2 line segments with common endpoints" << nl;
        segment_comparer_type cmp{{0, 0}};

        test_line_segments_are_equal(cmp, {1, 1}, {1, 0}, {1, 0}, {1, -1});
        test_line_segments_are_equal(cmp, {1, 1}, {1, 0}, {1, 0}, {1, 1});
        test_line_segment_is_closer(cmp, {2, 0}, {1, 1}, {2, 1}, {2, 0});
        test_line_segment_is_closer(cmp, {2, 1}, {2, 0}, {2, 0}, {3, 1});
    }

    void T2() {
        cout << "Compare angle with 2 points in general position" << nl;
        angle_comparer_type cmp{{0, 0}};

        REQUIRE(cmp({0, 1}, {1, 1}));
        REQUIRE_FALSE(cmp({1, 1}, {0, 1}));

        REQUIRE(cmp({1, 1}, {1, -1}));
        REQUIRE_FALSE(cmp({1, -1}, {1, 1}));

        REQUIRE(cmp({1, 0}, {-1, -1}));
        REQUIRE_FALSE(cmp({-1, -1}, {1, 0}));

        REQUIRE(cmp({0, 1}, {0, -1}));
        REQUIRE_FALSE(cmp({0, -1}, {0, 1}));
    }

    void T3() {
        cout << "Compare angle with 2 points if they are collinear with the origin" << nl;
        angle_comparer_type cmp{{0, 0}};

        REQUIRE(cmp({1, 0}, {2, 0}));
        REQUIRE_FALSE(cmp({2, 0}, {1, 0}));

        REQUIRE_FALSE(cmp({1, 0}, {1, 0}));
        REQUIRE_FALSE(cmp({0, 0}, {0, 0}));
    }

    void T4() {
        cout << "Calculate visibility polygon with no line segments" << nl;
        std::vector<segment_type> segments;
        auto poly = visibility_polygon(vector_type{0, 0}, {segments});
        REQUIRE(poly.p.empty());
    }

    void T5() {
        cout << "Calculate visibility polygon with no obstaces apart from the boundary" << nl;
        std::vector<segment_type> segments{
                {{-250, -250}, {-250, 250}},
                {{-250, 250},  {250,  250}},
                {{250,  250},  {250,  -250}},
                {{250,  -250}, {-250, -250}}
        };

        auto poly = visibility_polygon(vector_type{0, 0}, {segments});
        REQUIRE(poly.p.size() == 4);
        REQUIRE(equal(poly.p[0], {250, 250}));
        REQUIRE(equal(poly.p[1], {250, -250}));
        REQUIRE(equal(poly.p[2], {-250, -250}));
        REQUIRE(equal(poly.p[3], {-250, 250}));
    }

    void T6() {
        cout << "Calculate visibility polygon with a polyline as an obstacle" << nl;
        std::vector<segment_type> segments{
                {{-250, -250}, {-250, 250}},
                {{-250, 250},  {250,  250}},
                {{250,  250},  {250,  -250}},
                {{250,  -250}, {-250, -250}},

                {{-50,  50},   {50,   50}},
                {{50,   50},   {50,   -50}},
        };

        auto poly = visibility_polygon(vector_type{0, 0}, {segments});
        REQUIRE(poly.p.size() == 6);
        REQUIRE(equal(poly.p[0], {50, 50}));
        REQUIRE(equal(poly.p[1], {50, -50}));
        REQUIRE(equal(poly.p[2], {250, -250}));
        REQUIRE(equal(poly.p[3], {-250, -250}));
        REQUIRE(equal(poly.p[4], {-250, 250}));
        REQUIRE(equal(poly.p[5], {-50, 50}));
    }

    void T7() {
        cout << "Calculate visibility polygon with a convex polygon as an obstacle" << nl;
        std::vector<segment_type> segments{
                {{-250, -250}, {-250, 250}},
                {{-250, 250},  {250,  250}},
                {{250,  250},  {250,  -250}},
                {{250,  -250}, {-250, -250}},

                {{-50,  50},   {50,   50}},
                {{50,   50},   {50,   100}},
                {{50,   100},  {-50,  100}},
                {{-50,  100},  {-50,  50}},
        };

        auto poly = visibility_polygon(vector_type{0, 0}, {segments});
        REQUIRE(poly.p.size() == 6);
        REQUIRE(equal(poly.p[0], {50, 50}));
        REQUIRE(equal(poly.p[1], {250, 250}));
        REQUIRE(equal(poly.p[2], {250, -250}));
        REQUIRE(equal(poly.p[3], {-250, -250}));
        REQUIRE(equal(poly.p[4], {-250, 250}));
        REQUIRE(equal(poly.p[5], {-50, 50}));
    }

    void T8() {
        cout << "Calculate visibility polygon with a concave polygon as an obstacle" << nl;
        std::vector<segment_type> segments{
                {{-250, -250}, {-250, 250}},
                {{-250, 250},  {250,  250}},
                {{250,  250},  {250,  -250}},
                {{250,  -250}, {-250, -250}},

                {{-50,  50},   {0,    100}},
                {{0,    100},  {50,   50}},
                {{50,   50},   {0,    200}},
                {{0,    200},  {-50,  50}},
        };

        auto poly = visibility_polygon(vector_type{0, 0}, {segments});
        REQUIRE(poly.p.size() == 7);
        REQUIRE(equal(poly.p[0], {0, 100}));
        REQUIRE(equal(poly.p[1], {50, 50}));
        REQUIRE(equal(poly.p[2], {250, 250}));
        REQUIRE(equal(poly.p[3], {250, -250}));
        REQUIRE(equal(poly.p[4], {-250, -250}));
        REQUIRE(equal(poly.p[5], {-250, 250}));
        REQUIRE(equal(poly.p[6], {-50, 50}));
    }

    void T9() {
        // TODO
    }

    void run() {
        T0();
        T1();
        T2();
        T3();
        T4();
        T5();
        T6();
        T7();
        T8();
    }
}
// _____________________________________________________
// TEST CASES

void test() {
    TestCases::run();
}


void solve() {
    // TODO find an improved method for finding the visibility polygon
    int n;
    cin >> n;
    PolyP a;
    a.p.resize(n);
    for (auto& p: a.p) cin >> p.x >> p.y;
    const Poly& as = convert(a);
    point2d g, t;
    cin >> g.x >> g.y >> t.x >> t.y;
    PolyP vis = visibility_polygon(t, as);

    for (auto& p: vis.p) cout << p << ", ";
    cout << vis.p[0] << nl;
    const Poly& viss = convert(vis);

//    map<pii, double> m;
//    for (auto& p: a.p) {
//        const auto& v = PolyPDistanceVec(p, viss);
//        if (doIntersect(v, as)) m[{p.x, p.y}] = INF;
//        else m[{p.x, p.y}] = len(v);
//    }
//    for (auto& p: m) {
//        cout << "(" << p.F.F << ", " << p.F.S << ")" << ": " << p.S << nl;
//    }
}

int32_t main() {
    fast
    freopen("./in", "r", stdin);
//    test();
    solve();
}
