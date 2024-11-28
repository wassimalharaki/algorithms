/**
 * 23:39:22 9/18/24
 * D
 */
// ./ICPC/Geometry/Medium/D.cpp
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
const double eps = 1e-12l;

bool equal(ftype a, ftype b) {
    return abs(a - b) <= max(abs(a), abs(b)) * eps;
}

bool lessThan(ftype a, ftype b) {
    return (b - a) > max(abs(a), abs(b)) * eps;
}

struct P {
    ftype x, y;
    int h;
    P(): x{0}, y{0}, h{0} {}
    P(ftype x, ftype y, int h): x(x), y(y), h(h) {}
    bool operator==(const P& p) const {
        return equal(abs(), p.abs());
    }
    bool operator<(const P& p) const {
        return lessThan(abs(), p.abs());
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
    [[nodiscard]] ftype abs() const {
        return sqrt(norm());
    }
    [[nodiscard]] P normalize() const {
        ftype v = abs();
        if (equal(0, v)) return *this;
        return *this / v;
    }
    [[nodiscard]] P normal() const {
        return { -y, x, 0 };
    }
    [[nodiscard]] double arg() const {
        return atan2(y, x);
    }
};
P operator*(ftype a, P b) {
    return b * a;
}
istream& operator>>(istream& in, P& p) {
    int x, y, h;
    in >> x >> y >> h;
    p.x = x;
    p.y = y;
    p.h = h;
    return in;
}
ostream& operator<<(ostream& out, const P& p) {
    return out << "(" << p.x << ", " << p.y << ")";
}


void solve() {
    int n, x, y;
    cin >> n >> x >> y;
    P morty{(ftype)x, (ftype)y, 0};
    vec<P> a(n);
    for (auto& p: a) cin >> p;
    auto comp = [](double p0, double p1) -> bool {
        return lessThan(p0, p1);
    };
    map<double, vec<P>, decltype(comp)> m;
    for (auto& p: a) {
        P v = (p - morty);
        m[v.arg()].push_back(p);
    }
    int cnt = 0;
    for (auto& pm: m) {
        auto ps = pm.S;
        sort(all(ps), [&](const P& p0, const P& p1) {
            return (p0 - morty) < (p1 - morty);
        });
        int curh = 0;
        for (auto& p: ps) {
            if (curh < p.h) {
                curh = p.h;
                cnt++;
            }
        }
    }
    cout << cnt << nl;
}

int32_t main() {
    fast
    const string NAME{"jaguar"};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
    }
    LT solve();
}
