/**
 * 12:17:16 9/21/24
 * O
 */
// ./ICPC/Geometry/Medium/O.cpp
#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;
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
using ftype = double;
const double eps = 1e-9l, inf = 1e13l;

bool equal(double a, double b) {
    return abs(a - b) <= max(abs(a), abs(b)) * eps;
}

bool lessThan(double a, double b) {
    return (b - a) > max(abs(a), abs(b)) * eps;
}

struct P {
    int index;
    ftype x, y;
    P(): x{0}, y{0}, index{-1} {}
    P(ftype x, ftype y, int i): x(x), y(y), index(i) {}
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
    ftype operator*(const P& p) const { // Cross product
        return x * p.y - y * p.x;
    }
    P operator/(ftype t) const {
        return P(*this) /= t;
    }
    [[nodiscard]] ftype dot(const P& p) const {
        return x * p.x + y * p.y;
    }
    [[nodiscard]] ftype norm() const {
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
        return { -y, x, index };
    }
};
P operator*(ftype a, P b) {
    return b * a;
}
istream& operator>>(istream& in, P& p) {
    int x, y;
    in >> x >> y;
    p.x = x;
    p.y = y;
    return in;
}
ostream& operator<<(ostream& out, const P& p) {
    return out << "(" << p.x << ", " << p.y << ")";
}
using seg = pair<P, P>;

P segmentIntersection(const seg& s0, const seg& s1) {
    P p0 = s0.F, p1 = s0.S, p2 = s1.F, p3 = s1.S;
    double p0_x = p0.x, p0_y = p0.y;
    double p1_x = p1.x, p1_y = p1.y;
    double p2_x = p2.x, p2_y = p2.y;
    double p3_x = p3.x, p3_y = p3.y;

    double s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    double s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    double t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);
    return {p0_x + (t * s1_x), p0_y + (t * s1_y), -1};
}

double distance(const P& p1, const P& p2) {
    return hypot( p2.y-p1.y, p2.x-p1.x );
}

double heron(const P& p1, const P& p2, const P& p3, bool mx) {
    double area = inf;
    double a = distance( p1, p2 );
    double b = distance( p1, p3 );
    double c = distance( p3, p2 );
    double s = (a+b+c)/2.0;
    area = sqrt( s*(s-a)*(s-b)*(s-c) );
    // NOTE: modify this if you want to calculate for the maximum as well
    if (equal(area, 0)) {
        return mx ? -1 : inf;
    }
    return area;
}

void solve() {
    int n;
    cin >> n;
    vec<seg> tmp(n);
    for (int i = 0; i < n; i++) {
        int x, y;
        cin >> x >> y;
        tmp[i].F = P(x, y, -1);
        cin >> x >> y;
        tmp[i].S = P(x, y, -1);
    }
    auto cmp = [](const P& p0, const P& p1) -> bool {
        return p0 < p1;
    };
    set<P, decltype(cmp)> pnts;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            pnts.insert(segmentIntersection(tmp[i], tmp[j]));
        }
    }
    // Already sorted by increasing X
    vec<P> a(pnts.begin(), pnts.end());
    n = (int)a.size();
    for (int i = 0; i < n; i++) a[i].index = i;
    int k = 0;
    vec<seg> segs(n * (n - 1) / 2);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            segs[k++] = {a[i], a[j]};
        }
    }
    sort(segs.begin(), segs.end(), [&](const auto& s0, const auto& s1) {
        P p0 = s0.S - s0.F, p1 = s1.S - s1.F;
        double a0 = atan2(p0.y, p0.x), a1 = atan2(p1.y, p1.x);
        return lessThan(a0, a1);
    });
    double mn = inf, mx = 0;
    for (auto &s : segs) {
        double area;
        // Look at the most extreme points (0 and n-1)
        // to look for the largest triangle
//        if (s.F != a.back() and s.S != a.back()) {
//            area = heron(s.F, s.S, a.back(), true);
//            if (lessThan(mx, area)) mx = area;
//        }
//        if (s.F != a[0] and s.S != a[0]) {
//            area = heron(s.F, s.S, a[0], true);
//            if (lessThan(mx, area)) mx = area;
//        }

        int index1 = s.F.index, index2 = s.S.index;
        if (index1 > 0 and s.F != a[index1 - 1] and s.S != a[index1 - 1]) {
            area = heron(s.F, s.S, a[index1 - 1], false);
            if (lessThan(area, mn)) mn = area;
        }
        if (index1 < n - 1 and s.F != a[index1+1] and s.S != a[index1+1]) {
            area = heron(s.F, s.S, a[index1+1], false);
            if (lessThan(area, mn)) mn = area;
        }
        if (index2 > 0 and s.F != a[index2 - 1] and s.S != a[index2 - 1]) {
            area = heron(s.F, s.S, a[index2 - 1], false);
            if (lessThan(area, mn)) mn = area;
        }
        if (index2 < n - 1 and s.F != a[index2+1] and s.S != a[index2+1]) {
            area = heron(s.F, s.S, a[index2+1], false );
            if (lessThan(area, mn)) mn = area;
        }

        // As we rotate the calipers, these two points
        // swap positions in the array
        a[index1] = s.S;
        s.S.index = index1;
        a[index2] = s.F;
        s.F.index = index2;
    }
    if (equal(mn, inf)) {
        cout << setprecision(0) << -1 << nl;
    } else {
        cout << setprecision(6) << mn << nl;
    }
}

int32_t main() {
    fast
    const string NAME{"bird"};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
    }
    cout << fixed;
    LT solve();
}
