/**
 * 08:58:20 11/8/24
 * AnimeshDecidesToSettleDown
 */
// ./ICPC/Geometry/Medium/I.cpp
#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;
#define int long long
#define uint unsigned int
#define double long double
#define fast ios_base::sync_with_stdio(false);cin.tie(NULL);cout<<fixed<<setprecision(25);
#define nl '\n'
#define all(v) v.begin(), v.end()
#define allr(v) v.rbegin(), v.rend()
#define NO (cout << "NO" << nl);
#define YES (cout << "YES" << nl);
#define F first
#define S second
#define pb push_back
#define emb emplace_back
#define INF LONG_LONG_MAX
#define MOD 1000000007ll
#define EPS 1e-9l
#define PI 3.14159265358979323846264338327950288L
#define pii pair<int, int>
#define X real()
#define Y imag()
#define vec vector
#define LT int T; cin >> T; while (T--)
template<typename T>
using vec2d = vector<vector<T>>;
template<typename T>
using ordered_set = tree<T, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;
using indexed_set = tree<int, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;

namespace Geometry {
    using ftype = double;
    using atype = double;
    const double eps = 1e-9l;
    const double pi = acos(-1.0l);
    const double e = exp(1.0L);
    const double inf = INF;
    const int iinf = LONG_LONG_MAX;
    const int prec = 25;

    ftype sqr(ftype v) { return v * v; }

    void setPrec(ostream& out = cout) {
        out << fixed << setprecision(prec);
    }

    //bool eq(const ftype& a, const ftype& b) {
    //    return a == b;
    //}
    bool eq(const ftype& a, const ftype& b) {
        return abs(a - b) <= eps;
    }

    //bool ls(const ftype& a, const ftype& b) {
    //    return a < b;
    //}
    bool ls(const ftype& a, const ftype& b) {
        return (b - a) > eps;
    }

    //bool lse(const ftype& a, const ftype& b) {
    //    return a <= b;
    //}
    bool lse(const ftype& a, const ftype& b) {
        return ls(a, b) or eq(a, b);
    }

    //bool gt(const ftype& a, const ftype& b) {
    //    return a > b;
    //}
    bool gt(const ftype& a, const ftype& b) {
        return not ls(a, b) and not eq(a, b);
    }

    // Returns angle with respect to the x-axis.
    atype angle(atype x, atype y) {
        return atan2(y, x);
    }
    atype torad(atype a) { return a / 180 * pi; }
    atype todeg(atype a) { return a / pi * 180; }

    int sign(const ftype& v) {
        if (ls(0, v)) return 1;
        else if (ls(v, 0)) return -1;
        return 0;
    }

    // Orientation of 3 points forming an angle
    enum class Orientation {
        l = 1,
        r = -1,
        c = 0
    };

    // 2D point
    struct P {
        ftype x, y;
        static P polar(const ftype& r, const atype& a) {
            return {r * cos(a), r * sin(a)};
        }
        P(): x{0}, y{0} {}
        P(ftype x, ftype y): x{x}, y{y} {}
        explicit P(const complex<ftype>& p): P(p.X, p.Y) {}
        explicit P(const pair<ftype, ftype>& p): P(p.F, p.S) {}
        P(const P& p): P(p.x, p.y) {}
        explicit operator complex<ftype>() const {
            return {x, y};
        }
        explicit operator pair<ftype, ftype>() const {
            return {x, y};
        }
        ftype dist(const P& p) const {
            return (*this - p).abs();
        }
        atype r() const {
            return hypot(x, y);
        }
        atype a() const {
            return angle(x, y);
        }
        atype ua(const P& p) const {
            // undirected angle
            double v0 = fabs(a() - p.a()), v1 = fabs(p.ap() - ap());
            if (ls(v0, v1)) return v0;
            else return v1;
        }
        atype ap() const {
            atype res = angle(x, y);
            if (res < 0) res += 2 * pi;
            return res;
        }
        P operator-() const {
            return *this * -1;
        }
        P operator+() const {
            return *this;
        }
        P operator+=(const P& p) {
            x += p.x;
            y += p.y;
            return *this;
        }
        P operator-=(const P& p) {
            x -= p.x;
            y -= p.y;
            return *this;
        }
        P operator*=(const ftype& v) {
            x *= v;
            y *= v;
            return *this;
        }
        P operator/=(const ftype& v) {
            x /= v;
            y /= v;
            return *this;
        }
        P operator%=(const ftype& v) {
            x = fmod(x, v);
            y = fmod(y, v);
            return *this;
        }
        P operator^=(const ftype& an) {
            return *this = rotateccw(an);
        }
        P operator|=(const ftype& an) {
            return *this = rotatecw(an);
        }
        P operator|(const ftype& an) {
            P res = *this;
            res |= an;
            return res;
        }
        P operator^(const ftype& an) {
            P res = *this;
            res ^= an;
            return res;
        }
        P operator+(const P& p) const {
            P res = *this;
            res += p;
            return res;
        }
        P operator-(const P& p) const {
            P res = *this;
            res -= p;
            return res;
        }
        P operator*(const ftype& v) const {
            P res = *this;
            res *= v;
            return res;
        }
        P operator/(const ftype& v) const {
            P res = *this;
            res /= v;
            return res;
        }
        P operator%(const ftype& v) const {
            P res = *this;
            res %= v;
            return res;
        }
        bool operator==(const P& p) const {
            return eq(x, p.x) and eq(y, p.y);
        }
        bool operator<(const P& p) const {
            return eq(p.x, x) ? ls(y, p.y) : ls(x, p.x);
        }
        bool operator<=(const P& p) const {
            return *this < p or *this == p;
        }
        ftype cross(const P& b, const P& c) const {
            return (b - *this).cross(c - *this);
        }
        Orientation orientation(const P& pb, const P& pc) const {
            const P& pa = *this;
            ftype d = (pb - pa).cross(pc - pa);
            return static_cast<Orientation>(
                    ls(0, d) - ls(d, 0)
            );
        }
        ftype cross(const P& p) const {
            return x * p.y - y * p.x;
        }
        ftype dot(const P& p) const {
            return x * p.x + y * p.y;
        }
        ftype norm() const {
            return dot(*this);
        }
        atype abs() const {
            return sqrt((atype)norm());
        }
        P truncate(ftype v) const {
            // returns a vector with norm v and having same direction
            ftype k = abs();
            if (sign(k) == 0) return *this;
            v /= k;
            return {x * v, y * v};
        }
        P normalize() const {
            return truncate(1);
        }
        P normal_l() const {
            return {-y, x};
        }
        P normal_r() const {
            return {y, -x};
        }
        P rotateccw(const atype& angle, const P& ref = {0, 0}) const {
            P res;
            ftype xs = x - ref.x, ys = y - ref.y;
            res.x = ref.x + (xs * cos(angle) - ys * sin(angle));
            res.y = ref.y + (xs * sin(angle) + ys * cos(angle));
            return res;
        }
        P rotatecw(const atype& angle, const P& ref = {0, 0}) const {
            return rotateccw(-angle, ref);
        }
        friend P operator*(const ftype& v, const P& p) {
            return p * v;
        }
        friend istream& operator>>(istream& in, P& p) {
            int xv, yv;
            in >> xv >> yv;
            p = P(xv, yv);
            return in;
        }
        friend ostream& operator<<(ostream& out, const P& p) {
            setPrec(out);
            return out << "(" << p.x << ", " << p.y << ")";
        }
    };

    // basic comparison functions
    auto angleCmp = [](const P& p0, const P& p1) -> bool {
        return ls(p0.a(), p1.a());
    };
    auto radiusCmp = [](const P& p0, const P& p1) -> bool {
        return ls(p0.r(), p1.r());
    };
    auto stdCmp = [](const P& p0, const P& p1) -> bool {
        return p0 < p1;
    };
    auto xCmp = [](const P& p0, const P& p1) -> bool {
        return ls(p0.x, p1.x);
    };
    auto yCmp = [](const P& p0, const P& p1) -> bool {
        return ls(p0.y, p1.y);
    };

    // Hash function
    struct p_hash {
        static __int128 splitmix128(__int128 x) {
            // gr = 0x9e3779b97f4a7c15f39cc0605cedc835
            // c0 = 0xbf58476d1ce4e5b9a3f7b72c1e3c9e3b
            // c1 = 0x94d049bb133111ebb4093822299f31d7
            __int128 gr = 0x9e3779b97f4a7c15, c0 = 0xbf58476d1ce4e5b9, c1 = 0x94d049bb133111eb;
            gr <<= 64;
            c0 <<= 64;
            c1 <<= 64;
            gr += 0xf39cc0605cedc835;
            c0 += 0xa3f7b72c1e3c9e3b;
            c1 += 0xb4093822299f31d7;
            x += gr;
            x = (x ^ (x >> 62)) * c0;
            x = (x ^ (x >> 59)) * c1;
            return x ^ (x >> 63);
        }

        size_t operator()(const P& p) const {
            static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
            __int128 x = (((__int128)p.x) << 64) | ((__int128)p.y);
            return splitmix128(x + FIXED_RANDOM);
        }
    };
}
using namespace Geometry;

namespace Geometry {
    pair<int, int> reduce(ftype v0, ftype v1) {
        // log(n)
        int a = (int)v0, b = (int)v1;
        if (a == 0) return {0, 1};
        if (b == 0) return {iinf, 1};
        int g = gcd(a, b);
        return {a / g, b / g};
    }

    enum class InterT {
        no = 0,
        yes = 1,
        inf = 2
    };

    struct VectorInter {
        InterT t;
        P p;
        ftype pr0, pr1;
    };

    struct V {
        // P normal form is a(x - x0) + b(y - y0) = 0.
        // Where (a, b) is a normal to the line.
        P s, d;

        static V normalForm(ftype a, ftype b, ftype c) {
            // ax + by + c = 0
            c *= -1;
            if (eq(0, a) and eq(0, b)) {
                // Invalid line
                throw exception();
            } else if (eq(0, a)) {
                return {{0, c / b}, {1, c / b}};
            } else if (eq(0, b)) {
                return {{c / a, 0}, {c / a, 1}};
            } else {
                return {{0, c / b}, {1, (c - a) / b}};
            }
        }

        V(): V({0, 0}, {0, 0}) {}
        V(const ftype& a, const ftype& b): s{0, b}, d{1, a + b} {
            // y = ax + b
        }
        V(const ftype& y1, const ftype& x1, const ftype& m): s{x1, y1}, d{m + x1, y1 - 1} {
            // y - y1 = (x1 - x) / m
        }
        V(const P& s, const P& d): s{s}, d{d} {}
        V(const ftype& x0, const ftype& y0, const ftype& x1, const ftype& y1): V({x0, y0}, {x1, y1}) {}
        P r() const {
            return d - s;
        }
        array<ftype, 3> normalForm() const {
            // ax + by + c = 0
            return {d.y - s.y, s.x - d.x, d.cross(s)};
        }
        atype a() const {
            return r().a();
        }
        ftype abs() const {
            return r().abs();
        }
        ftype norm() const {
            return r().norm();
        }
        ftype dist(const P& p) const {
            if (eq(norm(), 0)) return p.dist(s);
            ftype t = (p - s).dot(r()) / norm();
            if (ls(t, 0)) t = 0;
            if (ls(1, t)) t = 1;
            return petP(t).dist(p);
        }
        ftype ldist(const V& v) const {
            VectorInter inter = intersection(v);
            if (inter.t != InterT::no) return 0;
            ftype vs = ldist(v.s), vd = ldist(v.d);
            ftype ss = v.ldist(s), sd = v.ldist(d);
            array<ftype, 4> tmp{vs, vd, ss, sd};
            sort(tmp.begin(), tmp.end(), ls);
            return tmp[0];
        }
        ftype dist(const V& v) const {
            VectorInter inter = intersection(v);
            if (inter.t != InterT::no) return 0;
            ftype vs = dist(v.s), vd = dist(v.d);
            ftype ss = v.dist(s), sd = v.dist(d);
            array<ftype, 4> tmp{vs, vd, ss, sd};
            sort(tmp.begin(), tmp.end(), ls);
            return tmp[0];
        }
        ftype rdist(const P& p) const {
            if (eq(norm(), 0)) return p.dist(s);
            ftype t = (p - s).dot(r()) / norm();
            if (ls(t, 0)) t = 0;
            return petP(t).dist(p);
        }
        ftype ldist(const P& p) const {
            return project(p).dist(p);
        }
        P project(const P& p) const {
            // Perpendicular projection of the point on the vector
            const ftype& t = (p - s).dot(r()) / norm();
            return petP(t);
        }
        Orientation orientation(const P& p) const {
            return s.orientation(d, p);
        }
        P petX(const ftype& x) const {
            // NOTE THAT IF d.x == s.x, then any point with p.x = x is valid
            // Returns point on the line formed by the vector based on given x
            // x = s.x + t * (d.x - s.x)
            // (x - s.x) / (d.x - s.x) = t
            if (eq(d.x, s.x)) return {x, d.y};
            return petP((x - s.x) / (d.x - s.x));
        }
        P petY(const ftype& y) const {
            // NOTE THAT IF d.y == s.y, then any point with p.y = y is valid
            // Returns point on the line formed by the vector based on given y
            // y = s.y + t * (d.y - s.y)
            // (y - s.y) / (d.y - s.y) = t
            if (eq(d.y, s.y)) return {d.x, y};
            return petP((y - s.y) / (d.y - s.y));
        }
        P petP(const ftype& p) const {
            // Returns point on the line formed by the vector based on given percentage of size
            return s + p * (d - s);
        }
        P reflect(const P& p) const {
            const P& l1 = s, l2 = d;
            P z = p - l1, w = l2 - l1;
            complex<ftype> zc = {z.x, z.y},
                    wc = {w.x, w.y},
                    l1c = {l1.x, l1.y};
            complex<ftype> res = conj(zc / wc) * wc + l1c;
            return {res.X, res.Y};
        }
        VectorInter intersection(const V& v) const {
            P p0 = s, p1 = d, p2 = v.s, p3 = v.d;
            ftype p0_x = p0.x, p0_y = p0.y;
            ftype p1_x = p1.x, p1_y = p1.y;
            ftype p2_x = p2.x, p2_y = p2.y;
            ftype p3_x = p3.x, p3_y = p3.y;

            ftype s1_x, s1_y, s2_x, s2_y;
            s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
            s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

            ftype t0 = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
            ftype t1 = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);
            P p{p0_x + (t1 * s1_x), p0_y + (t1 * s1_y)};

            // Parallel
            if (isinf(t0) or isinf(t1)) {
                return {
                        InterT::no,
                        {inf, inf},
                        t0,
                        t1
                };
            }

            // Collinear
            if (isnan(t0) or isnan(t1)) {
                bool os = on(v.s), od = on(v.d);
                if (os or od) {
                    p = (os ? v.s : v.d);
                    return {
                            InterT::inf,
                            p,
                            getP(p).S,
                            (os ? 0.0l : 1.0l)
                    };
                }

                return {
                        InterT::no,
                        {inf, inf},
                        t0,
                        t1
                };
            }

            // Regular
            if (lse(0, t0) and lse(t0, 1) and lse(0, t1) and lse(t1, 1)) {
                // Collision detected
                return {
                        InterT::yes,
                        p,
                        t0,
                        t1
                };
            }

            return {
                    InterT::no,
                    p,
                    t0,
                    t1
            };
        }
        ftype slope() const {
            const P& tmp = r();
            return tmp.y / tmp.x;
        }
        pair<int, int> fslope() const {
            return reduce(d.y - s.y, d.x - s.x);
        }
        atype angle(const V& v) const {
            // Returns NAN if the lines are parallel or collinear
            return acos(dot(v) / (abs() * v.abs()));
        }
        pair<char, ftype> getP(const P& p) const {
            if (not onPath(p)) return {false, -1};
            // x = s.x + p * (d.x - s.x)
            // (x - s.x) / (d.x - s.x) = p
            // (y - s.y) / (d.y - s.y) = p
            if (eq(d.x, s.x)) {
                return {true, (p.y - s.y) / (d.y - s.y)};
            } else if (not eq(d.y, s.y)) {
                return {true, (p.x - s.x) / (d.x - s.x)};
            } else {
                return {p == s, 0};
            }
        }
        ftype left() const {
            return ls(d.x, s.x) ? d.x : s.x;
        }
        ftype right() const {
            return ls(d.x, s.x) ? s.x : d.x;
        }
        ftype bottom() const {
            return ls(d.y, s.y) ? d.y : s.y;
        }
        ftype top() const {
            return ls(d.y, s.y) ? s.y : d.y;
        }
        bool onPath(const P& p) const {
            if (eq(norm(), 0)) return p == s;
            return eq(r().cross(p - s), 0);
        }
        bool on(const P& p) const {
            // Not collinear or behind
            return onPath(p) and lse(left(), p.x) and lse(p.x, right());
        }
        bool onRay(const P& p) const {
            if (not onPath(p)) return false;
            if (lse(s.x, d.x)) return lse(s.x, p.x);
            return lse(p.x, s.x);
        }
        bool parallel(const V& v) const {
            return isnan(angle(v));
        }
        bool orthogonal(const V& v) const {
            return eq(dot(v), 0);
        }
        V operator-() const {
            return *this * -1;
        }
        V operator+() const {
            return *this;
        }
        V operator+=(const V& v) {
            s += v.s;
            d += v.d;
            return *this;
        }
        V operator-=(const V& v) {
            s -= v.s;
            d -= v.d;
            return *this;
        }
        V operator*=(const ftype& v) {
            s *= v;
            d *= v;
            return *this;
        }
        V operator/=(const ftype& v) {
            s /= v;
            d /= v;
            return *this;
        }
        V operator%=(const ftype& v) {
            s %= v;
            d %= v;
            return *this;
        }
        V operator^=(const ftype& an) {
            return *this = rotateccw(an);
        }
        V operator|=(const ftype& an) {
            return *this = rotatecw(an);
        }
        V operator|(const ftype& an) {
            V res = *this;
            res |= an;
            return res;
        }
        V operator^(const ftype& an) {
            V res = *this;
            res ^= an;
            return res;
        }
        V operator+(const V& p) const {
            V res = *this;
            res += p;
            return res;
        }
        V operator-(const V& p) const {
            V res = *this;
            res -= p;
            return res;
        }
        V operator*(const ftype& v) const {
            V res = *this;
            res *= v;
            return res;
        }
        V operator/(const ftype& v) const {
            V res = *this;
            res /= v;
            return res;
        }
        V operator%(const ftype& v) const {
            V res = *this;
            res %= v;
            return res;
        }
        bool operator==(const V& v) const {
            return s == v.s and d == v.d;
        }
        ftype cross(const V& v) const {
            return r().cross(v.r());
        }
        ftype dot(const V& v) const {
            return r().dot(v.r());
        }
        V truncate(ftype v) const {
            // returns a vector with norm v and having same direction
            ftype k = abs();
            if (sign(k) == 0) return *this;
            v /= k;
            return {s * v, d * v};
        }
        V normalize() const {
            return truncate(1);
        }
        V ndrotateccw(const atype& angle, const P& ref = {0, 0}) {
            V res{*this};
            res.d = res.d.rotateccw(angle, ref);
            res.s = res.s.rotateccw(angle, ref);
            return res;
        }
        V rotateccw(const atype& angle) const {
            V res{*this};
            res.d = res.d.rotateccw(angle, res.s);
            return res;
        }
        V rotatecw(const atype& angle) const {
            return rotateccw(-angle);
        }
        friend V operator*(const ftype& val, const V& v) {
            return v * val;
        }
        friend istream& operator>>(istream& in, V& v) {
            return in >> v.s >> v.d;
        }
        friend ostream& operator<<(ostream& out, const V& v) {
            return out << v.s << " -> " << v.d;
        }
    };
}

// Basic half-plane struct.
struct Halfplane {
    // 'p' is a passing point of the line and 'pq' is the direction vector of the line.
    P p, pq;
    double angle;

    Halfplane() {}
    Halfplane(const P& a, const P& b) : p(a), pq(b - a) {
        angle = atan2l(pq.y, pq.x);
    }

    // Check if point 'r' is outside this half-plane.
    // Every half-plane allows the region to the LEFT of its line.
    bool out(const P& r) const {
        return ls(pq.cross(r - p), 0);
    }

    // Comparator for sorting.
    bool operator < (const Halfplane& s) const {
        return ls(angle, s.angle);
    }

    // Intersection point of the lines of two half-planes. It is assumed they're never parallel.
    friend P inter(const Halfplane& s, const Halfplane& t) {
        double alpha = (t.p - s.p).cross(t.pq) / s.pq.cross(t.pq);
//        assert(not isnan(alpha));
        return s.p + (s.pq * alpha);
    }
};

// might contain duplicate points
// nlog(n)
// Actual algorithm
vector<P> hp_intersect(vector<Halfplane>& H) {

    P box[4] = {  // Bounding box in CCW order
            P(inf, inf),
            P(-inf, inf),
            P(-inf, -inf),
            P(inf, -inf)
    };

    for(int i = 0; i<4; i++) { // Add bounding box half-planes.
        Halfplane aux(box[i], box[(i+1) % 4]);
        H.push_back(aux);
    }

    // Sort by angle and start algorithm
    sort(H.begin(), H.end());
    deque<Halfplane> dq;
    int len = 0;
    for(int i = 0; i < (int)(H.size()); i++) {

        // Remove from the back of the deque while last half-plane is redundant
        while (len > 1 && H[i].out(inter(dq[len-1], dq[len-2]))) {
            dq.pop_back();
            --len;
        }

        // Remove from the front of the deque while first half-plane is redundant
        while (len > 1 && H[i].out(inter(dq[0], dq[1]))) {
            dq.pop_front();
            --len;
        }

        // Special case check: Parallel half-planes
        if (len > 0 && eq(H[i].pq.cross(dq[len-1].pq), 0)) {
            // Opposite parallel half-planes that ended up checked against each other.
            if (H[i].pq.dot(dq[len-1].pq) < 0.0)
                return vector<P>();

            // Same direction half-plane: keep only the leftmost half-plane.
            if (H[i].out(dq[len-1].p)) {
                dq.pop_back();
                --len;
            }
            else continue;
        }

        // Add new half-plane
        dq.push_back(H[i]);
        ++len;
    }

    // Final cleanup: Check half-planes at the front against the back and vice-versa
    while (len > 2 && dq[0].out(inter(dq[len-1], dq[len-2]))) {
        dq.pop_back();
        --len;
    }

    while (len > 2 && dq[len-1].out(inter(dq[0], dq[1]))) {
        dq.pop_front();
        --len;
    }

    // Report empty intersection if necessary
    if (len < 3) return vector<P>();

    // Reconstruct the convex polygon from the remaining half-planes.
    vector<P> ret(len);
    for(int i = 0; i+1 < len; i++) {
        ret[i] = inter(dq[i], dq[i+1]);
    }
    ret.back() = inter(dq[len-1], dq[0]);
    return ret;
}

double area(const vec<P>& v) {
    double res = 0;
    for (int i = 0; i < v.size(); i++) {
        res += v[i].cross(v[(i + 1) % v.size()]);
    }
    return abs(res) / 2;
}

void solve() {
    int n;
    cin >> n;
    vec<P> ps(n);
    for (P& p: ps) cin >> p;

    ftype total = area(ps);
    vec<Halfplane> h;
    for (int i = 0; i < n; i++) h.emplace_back(ps[i], ps[(i + 1) % n]);
    ftype all = area(hp_intersect(h));

    auto ans = all / total;
    if (eq(total, 0)) {
        cout << 0 << nl;
    } else {
        cout << ans << nl;
    }
}

int32_t main() {
    fast
    const string NAME{" "};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
        freopen((NAME + ".out").c_str(), "w", stdout);
    }
    solve();
}
