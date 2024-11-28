/**
 * 18:51:09 9/11/24
 * P
 */
// ./ICPC/Geometry/Medium/P.cpp
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

bool equal(double a, double b) {
    return abs(a - b) <= max(abs(a), abs(b)) * EPS;
}

bool lessThan(double a, double b) {
    return (b - a) > max(abs(a), abs(b)) * EPS;
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
        if (::abs(v) < EPS) return *this;
        return *this / v;
    }
    [[nodiscard]] P normal() const {
        return { -y, x };
    }
};
P operator*(ftype a, P b) {
    return b * a;
}

double signed_area_parallelogram(P& p1, P& p2, P& p3) {
    return (p2 - p1) * (p3 - p2);
}

double triangle_area(P& p1, P& p2, P& p3) {
    return abs(signed_area_parallelogram(p1, p2, p3)) / 2.0;
}

bool collinear(P& p0, P& p1, P& p2) {
    return equal(0, triangle_area(p0, p1, p2));
}

void solve() {
    int n;
    cin >> n;
    vec<P> p(n);
    for (auto& i: p) {
        int x, y;
        cin >> x >> y;
        i.x = x, i.y = y;
    }
    int cnt = n * (n - 1) * (n - 2) * (n - 3) / 24; // nC4
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            int collinear_cnt = 0, future_collinear_cnt = 0;
            for (int k = 0; k < n; k++) {
                if (k == i or k == j) continue;
                if (collinear(p[i], p[j], p[k])) {
                    collinear_cnt++;
                    if (k > j) future_collinear_cnt++;
                }
            }

            // Removes the cases where a collinear point (after j) forms a degenerate quadrilateral with non-collinear points.
            cnt -= future_collinear_cnt * (n - 2 - collinear_cnt);

            // Removes fully degenerate quadrilaterals where all four points are collinear.
            cnt -= future_collinear_cnt * (future_collinear_cnt - 1) / 2;
        }
    }
    cout << cnt << nl;
}

int32_t main() {
    fast
//    const string NAME{"zoser"};
//    if (NAME != " ") {
//        freopen((NAME + ".in").c_str(), "r", stdin);
//    }
    solve();
}
