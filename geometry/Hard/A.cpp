/**
 * 14:44:59 12/4/24
 * A
 */
// ./algorithms/geometry/Hard/A.cpp
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
#define INF LONG_LONG_MAX
#define pii pair<int, int>
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
//    const double e = exp(1.0L);
    const double inf = 1e18l;
    const int mod = 1e9 + 7;
    const int iinf = LONG_LONG_MAX;
    int32_t prec = 25;
    bool useDoubleIn = false;

    ftype sqr(ftype v) { return v * v; }

    void setPrec(ostream& out = cout) {
        out << fixed << setprecision(prec);
    }

    //bool eq(const ftype& a, const ftype& b) {
    //    return a == b;
    //}
    bool eq(const ftype& a, const ftype& b) {
        return abs(a - b) < eps;
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

    //bool bt(const ftype& a, const ftype& b, const ftype& c) {
    //    return a < b and b < c;
    //}
    bool bt(const ftype& a, const ftype& b, const ftype& c) {
        return ls(a, b) and ls(b, c);
    }

    //bool bte(const ftype& a, const ftype& b, const ftype& c) {
    //    return a <= b and b <= c;
    //}
    bool bte(const ftype& a, const ftype& b, const ftype& c) {
        return lse(a, b) and lse(b, c);
    }

    ftype max(const ftype& a, const ftype& b) {
        return ls(a, b) ? b : a;
    }

    int max(int a, int b) {
        return a < b ? b : a;
    }

    ftype min(const ftype& a, const ftype& b) {
        return ls(a, b) ? a : b;
    }

    int min(int a, int b) {
        return a < b ? a : b;
    }

    optional<array<ftype, 2>> quadratic(ftype a, ftype b, ftype c) {
        ftype discriminant = b * b - 4 * a * c;
        if (ls(discriminant, 0)) {
            return nullopt;
        }
        ftype sqrtD = sqrt(discriminant);
        return array<ftype, 2>{(-b + sqrtD) / (2 * a), (-b - sqrtD) / (2 * a)};
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

    pair<int, int> reduce(ftype v0, ftype v1) {
        // log(n)
        int a = (int)v0, b = (int)v1;
        if (a == 0) return {0, sign(b)};
        if (b == 0) return {sign(a) * iinf, 0};
        int g = gcd(a, b);
        return {a / g, b / g};
    }

    pair<int, int> reducef(ftype v0, ftype v1) {
        // log(n)
        return reduce(v0 * 1e6l, v1 * 1e6l);
    }

    // 2D point
    struct P {
        ftype x, y;
        static P polar(const ftype& r, const atype& a) {
            return {r * cos(a), r * sin(a)};
        }
        P(): x{0}, y{0} {}
        P(ftype x, ftype y): x{x}, y{y} {}
        explicit P(const complex<ftype>& p): P(p.real(), p.imag()) {}
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
        ftype dist2(const P& p) const {
            return (*this - p).norm();
        }
        atype r() const {
            return hypot(x, y);
        }
        atype a() const {
            return Geometry::angle(x, y);
        }
        atype ua(const P& p) const {
            // undirected angle
            double v0 = fabs(a() - p.a()), v1 = fabs(p.ap() - ap());
            if (ls(v0, v1)) return v0;
            else return v1;
        }
        atype ap() const {
            atype res = Geometry::angle(x, y);
            if (res < 0) res += 2 * pi;
            return res;
        }
        atype angle(const P& p) const {
            return acos(dot(p) / (abs() * p.abs()));
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
        bool ycmp(const P& p) const {
            return eq(p.y, y) ? ls(x, p.x) : ls(y, p.y);
        }
        bool operator<=(const P& p) const {
            return *this < p or *this == p;
        }
        ftype slope() const {
            return y / x;
        }
        pair<int, int> fslope() const {
            return reduce(y, x);
        }
        pair<int, int> inorm() const {
            return reduce(x, y);
        }
        ftype cross(const P& b, const P& c) const {
            return (b - *this).cross(c - *this);
        }
        bool is_point_in_angle(P b, const P& a, P c) const {
            // does point p lie in angle <bac
            assert(c.orientation(a, b) != 0);
            if (a.orientation(c, b) < 0) swap(b, c);
            return a.orientation(c, *this) >= 0 && a.orientation(b, *this) <= 0;
        }
        // 1 this is to right of ab, 0 collinear, -1 to left
        int orientation(const P& a, const P& b) const {
            // orientation of this with respect to AB
            const P& c = *this;
            ftype d = (b - c).cross(a - c);
            return sign(d);
        }
        ftype cross(const P& p) const {
            return x * p.y - y * p.x;
        }
        ftype dot(const P& p) const {
            return x * p.x + y * p.y;
        }
        ftype dot(const P& a, const P& b) const {
            return (a - *this).dot(b - *this);
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
            return *this * v;
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
        P normal_lfs() const {
            return P(fslope()).normal_l();
        }
        P normal_rfs() const {
            return P(fslope()).normal_r();
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
            if (useDoubleIn) {
                return in >> p.x >> p.y;
            }

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

    P inf_pt{inf, inf};

    double heron(const P& p1, const P& p2, const P& p3) {
        double a = p1.dist(p2);
        double b = p1.dist(p3);
        double c = p3.dist(p2);
        double s = (a+b+c)/2.0;
        double area = sqrt( s*(s-a)*(s-b)*(s-c) );
        return area;
    }

    bool isCollinear(const P& a, const P& b, const P& c) {
        return eq(0, heron(a, b, c));
    }
}
using namespace Geometry;
namespace Geometry {
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
        // Point normal form is a(x - x0) + b(y - y0) = 0.
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
        bool shares(const V& v) const {
            return s == v.s or s == v.d or d == v.s or d == v.d;
        }
        int sameDir(const V& v) const { // -1 opposite, 0 collinear, 1 same
            return sign(dot(v));
        }
        bool vertical() const {
            return eq(s.x, d.x);
        }
        V push(const V& v) const {
            return {s + v.r(), d + v.r()};
        }
        P r() const {
            return d - s;
        }
        P midpoint() const {
            return 0.5l * (s + d);
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
        ftype mx_dist(const P& p) const {
            return max(p.dist(s), p.dist(d));
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
        int orientation(const P& p) const {
            return p.orientation(s, d);
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
            return {res.real(), res.imag()};
        }

        V normal_l() const {
            P t = r().normal_l();
            return {s, s + t};
        }

        V normal_r() const {
            P t = r().normal_r();
            return {s, s + t};
        }

        void assign(const V& v) {
            s = v.d;
            d = v.d;
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
                if (lp() == v.rp()) {
                    p = lp();
                    return {
                            InterT::yes,
                            p,
                            getP(p).S,
                            v.getP(p).S
                    };
                }
                if (rp() == v.lp()) {
                    p = rp();
                    return {
                            InterT::yes,
                            p,
                            getP(p).S,
                            v.getP(p).S
                    };
                }
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

            // Regular, modify this for rays or lines
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
            return r().slope();
        }
        pair<int, int> fslope() const {
            return r().fslope();
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
        P& lp() {
            return s < d ? s : d;
        }
        P& rp() {
            return s < d ? d : s;
        }
        const P& lp() const {
            return s < d ? s : d;
        }
        const P& rp() const {
            return s < d ? d : s;
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
            return isCollinear(s, d, p);
        }
        bool on(const P& p) const {
            // Not collinear or behind
            if (not onPath(p)) return false;
            if (vertical()) return lse(bottom(), p.y) and lse(p.y, top());
            return lse(left(), p.x) and lse(p.x, right());
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
        bool lslope(const V& v) const {
            pii a = fslope(), b = v.fslope();
            return a.F * b.S < a.S * b.F;
        }
        bool operator==(const V& v) const {
            return s == v.s and d == v.d;
        }
        bool operator<(const V& v) const {
            if (lp() == v.lp()) {
                return rp() < v.rp();
            }
            return lp() < v.lp();
        }
        // To the left of the vector
        bool operator<(const P& p) const {
            return ls(0, r().cross(p - s));
        }
        // On the vector
        bool operator==(const P& p) const {
            return on(p);
        }
        // To the right of the vector
        bool operator>(const P& p) const {
            return ls(r().cross(p - s), 0);
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
            return out << v.s << ", " << v.d;
        }
        bool collinear(const V& v) const {
            return isCollinear(s, d, v.s) and isCollinear(s, d, v.d);
        }
    };
}

int N;
void dfs(vec2d<char>& adj, vec<char>& marked, int n, int vert, int start, int& count) {
    marked[vert] = true;
    if (n == 0) {
        marked[vert] = false;
        if (adj[vert][start] and adj[start][vert]) {
            count++;
        }
        return;
    }
    for (int i = 0; i < N; i++) {
        if (not marked[i] and adj[vert][i]) {
            dfs(adj, marked, n - 1, i, start, count);
        }
        marked[vert] = false;
    }
}

int C3(int n) { return n * (n - 1) * (n - 2) / 6; }

int countCycles(vec2d<char>& adj, int n) {
    vec<char> marked(N);
    int count = 0;
    for (int i = 0; i <= N - n; i++) {
        dfs(adj, marked, n - 1, i, i, count);
        marked[i] = true;
    }
    return count >> 1;
}

void solve() {
    int n;
    cin >> n;
    vec<V> a(n);
    for (V& v: a) cin >> v;
    map<P, set<int>> m;
    auto countZeroCycles = [&]() -> int {
        int cnt = 0;
        for (auto& [_, s]: m) {
            int k = (int)s.size();
            cnt += C3(k);
        }
        return cnt;
    };

    vec2d<char> adj(n, vec<char>(n));
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (i == j) continue;
            auto inter = a[i].intersection(a[j]);
            if (inter.t == InterT::yes) {
                adj[i][j] = adj[j][i] = true;
                m[inter.p].insert(i);
                m[inter.p].insert(j);
            }
        }
    }
    N = n;
    int res = countCycles(adj, 3) - countZeroCycles();
    assert(res >= 0);
    cout << res << nl;
}

int32_t main() {
    fast
    const string NAME{" "};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
        freopen((NAME + ".out").c_str(), "w", stdout);
    }
    LT solve();
}
