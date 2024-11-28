/**
 * 12:32:17 9/15/24
 * J
 */
// ./ICPC/Geometry/Medium/J.cpp
#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;
#define int long long
#define uint unsigned int
#define double long double
#define fast ios_base::sync_with_stdio(false);cin.tie(NULL);cout<<setprecision(16);
#define nl '\n'
#define all(v) v.begin(), v.end()
#define NO (cout << "NO" << nl);
#define YES (cout << "YES" << nl);
#define F first
#define S second
#define INF LONG_LONG_MAX
#define MOD 1000000007ll
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

using ftype = double;
const double pi =  3.14159265358979323846264338327950288l;
const double eps = 1e-9l;

bool equal(double a, double b) {
    return abs(a - b) <= max(abs(a), abs(b)) * eps;
}

bool lessThan(double a, double b) {
    return (b - a) > eps;
}

struct P {
    ftype x, y;
    P(): x{0}, y{0} {}
    P(ftype x, ftype y): x(x), y(y) {}
    bool operator==(const P& p) const {
        return equal(x, p.x) and equal(y, p.y);
    }
    bool operator<(const P& p) const {
        return equal(x, p.x) ? lessThan(y, p.y) : lessThan(x, p.x);
    }
    P& operator+=(const P &t) {
        x += t.x;
        y += t.y;
        return *this;
    }
    P& operator-=(const P &t) {
        x -= t.x;
        y -= t.y;
        return *this;
    }
    P& operator*=(ftype t) {
        x *= t;
        y *= t;
        return *this;
    }
    P& operator/=(ftype t) {
        x /= t;
        y /= t;
        return *this;
    }
    P operator+(const P &t) const {
        return P(*this) += t;
    }
    P operator-(const P &t) const {
        return P(*this) -= t;
    }
    P operator*(ftype t) const {
        return P(*this) *= t;
    }
    double operator*(const P& p) const { // Cross product
        return x * p.y - y * p.x;
    }
    P operator/(ftype t) const {
        return P(*this) /= t;
    }
    [[nodiscard]] double dot(const P& p) const {
        return x * p.x + y * p.y;
    }
    [[nodiscard]] double norm() const {
        return dot(*this);
    }
    [[nodiscard]] double abs() const {
        return sqrt(norm());
    }
    [[nodiscard]] P normalize() const {
        double v = abs();
        if (::abs(v) < eps) return *this;
        return *this / v;
    }
    [[nodiscard]] P normal() const {
        return { -y, x };
    }
};
istream& operator>>(istream& in, P& p) {
    return in >> p.x >> p.y;
}
ostream& operator<<(ostream& out, const P& p) {
    return out << "(" << p.x << ", " << p.y << ")";
}
using seg = pair<P, P>;

// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines
// intersect the intersection point may be stored in second element of the pair.
pair<char, P> segmentIntersection(const seg& s0, const seg& s1) {
    P p0 = s0.F, p1 = s0.S, p2 = s1.F, p3 = s1.S;
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

    if (not lessThan(s, 0)
        and not lessThan(1, s)
        and not lessThan(t, 0)
        and not lessThan(1, t)) {
        // Collision detected
        return {1, {p0_x + (t * s1_x), p0_y + (t * s1_y)}};
    }

    return {0, {}}; // No collision
}

P reflect(const P& p, const seg& s) {
    const P& l1 = s.F, l2 = s.S;
    P z = p - l1, w = l2 - l1;
    complex<ftype> zc = {z.x, z.y},
                   wc = {w.x, w.y},
                   l1c = {l1.x, l1.y};
    complex<ftype> res = conj(zc / wc) * wc + l1c;
    return {res.X, res.Y};
}

bool doIntersect(int ref, const seg& s, const vec<seg>& ss) {
    for (int j = 0; j < (int)ss.size(); j++) {
        if (j != ref and segmentIntersection(ss[j], s).F) {
            return true;
        }
    }
    return false;
}

bool can_reach(int si, const vec<seg>& ss, const P& hole, const P& hole_ref) {
    P origin = {0, 0};
    seg s = ss[si], refp = {origin, hole_ref};
    auto [inter, inter_p] = segmentIntersection(s, refp);
    if (not inter) return false;
    return not doIntersect(si, {origin, inter_p}, ss)
       and not doIntersect(si, {inter_p, hole}, ss);
}

void solve() {
    int n, x, y;
    cin >> n >> x >> y;
    vec<seg> ss(n);
    for (auto& p: ss) cin >> p.F >> p.S;
    P hole, start;
    hole.x = x, hole.y = y;
    start.x = 0, start.y = 0;
    int cnt = 0;
    for (int i = 0; i < n; i++) {
        const seg& s = ss[i];
        if (can_reach(i, ss, hole, reflect(hole, s))) cnt++;
    }
    cout << cnt << nl;
}

int32_t main() {
    fast
    const string NAME{" "};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
//        freopen((NAME + ".out").c_str(), "w", stdout);
    }
    solve();
}
