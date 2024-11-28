/**
 * 12:27:31 10/6/24
 * V
 */
// ./ICPC/Geometry/OrganizedTemplates/V.cpp
#include "./point.cpp"

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
        P lp() const {
            return s < d ? s : d;
        }
        P rp() const {
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
        bool lslope(const V& v) const {
            pii a = fslope(), b = v.fslope();
            return a.F * b.S < a.S * b.F;
        }
        bool operator==(const V& v) const {
            return s == v.s and d == v.d;
        }
        bool operator<(const V& v) const {
            return s < v.s or s == v.s and d < v.d;
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
    };
}

//void solve() {
//    V v0{0, 0, 1, 1}, v1{2, 2, 3, 4};
//    cout << v1 << nl << v1.normal_l() << nl << v1.normal_r() << nl;
////    V vp0{0, 0, 1, 1}, vp1{0, -1, 1, 0};
////    V ve0{0, 0, 1, 1};
////    P p{0, 0}, p1{1, 1}, p2{-2, -2}, p3{-2, -3};
////    cout << v0.on(p) << ' ' << v0.onPath(p) << nl;
////    cout << v0.on(p1) << ' ' << v0.onPath(p1) << nl;
////    cout << v0.on(p2) << ' ' << v0.onPath(p2) << nl;
////    cout << v0.on(p3) << ' ' << v0.onPath(p3) << nl;
////    const auto& i = v0.intersection(v1);
////    cout << i.p << ' ' << (int)i.t << nl;
////    cout << vp0.intersection(vp1).p << nl;
////    cout << ve0.intersection(ve0).p << nl;
//}

//int32_t main() {
//    fast
//    solve();
//}
