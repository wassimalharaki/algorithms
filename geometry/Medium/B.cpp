/**
 * 15:15:28 9/27/24
 * Point
 */
// ./ICPC/Geometry/OrganizedTemplates/Point.cpp
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
    int32_t prec = 25;

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

    //bool gt(const ftype& a, const ftype& b) {
    //    return a > b;
    //}
    bool gt(const ftype& a, const ftype& b) {
        return not ls(a, b) and not eq(a, b);
    }

    ftype max(const ftype& a, const ftype& b) {
        return ls(a, b) ? b : a;
    }

    ftype min(const ftype& a, const ftype& b) {
        return ls(a, b) ? a : b;
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

using Node = pair<double, pii>;

double minDisToX(const P& from, const P& to, double X, double alpha, double beta) {
    if (eq(to.x, X)) return 0;

    P Vx(X - to.x, 0);
    double theta = Vx.angle(from - to);

    if (ls(theta, alpha)) {
        if (lse(0, Vx.cross(from - to))) Vx ^= alpha - theta;
        else Vx ^= theta - alpha;

    } else if (ls(beta, theta)) {
        if (lse(0, Vx.cross(from - to))) Vx ^= theta - beta;
        else Vx ^= beta - theta;
    }
    if (ls(abs(Vx.x), 0)) return inf;
    double m = (X - to.x) / Vx.x;
    if (ls(m, 0)) return inf;
    return (Vx * m).abs();
}

void solve() {
    int n, l, s, alpha_i, beta_i;
    cin >> n >> l >> s >> alpha_i >> beta_i;
    double alpha = torad(alpha_i), beta = torad(beta_i);
    vec<P> ps(n);
    vec2d<pair<int, double>> adj(n);
    for (P& p: ps) cin >> p;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double d = ps[i].dist(ps[j]);
            if (ls(d, l)) {
                adj[i].emplace_back(j, d);
                adj[j].emplace_back(i, d);
            }
        }
    }
    vec2d<double> dis(n + 1, vec<double>(n + 1));
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            dis[i][j] = inf;
        }
    }
    auto cmp = [](const Node& a, const Node& b) -> bool {
        if (eq(a.F, b.F)) return a.S < b.S;
        return ls(a.F, b.F);
    };
    priority_queue<Node, vec<Node>, decltype(cmp)> q;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double d2 = ps[i].dist(ps[j]);
            if (i == j or ls(l, d2) or ls(l, ps[i].x)) continue;
            double d = minDisToX(ps[j], ps[i], 0, alpha, beta);
            if (ls(l, d)) continue;
            dis[i][j] = d + d2;
            q.push({-(d + d2), {i, j}});
        }
    }
    double ans = inf;
    while (not q.empty()) {
        double cost = -q.top().F;
        auto [x, y] = q.top().S;
        q.pop();
        if (ls(ans, cost)) continue;
        double d = minDisToX(ps[x], ps[y], s, alpha, beta);

        if (ls(d, l)) {
            ans = min(ans, cost + d);
        }
        for (auto& p: adj[y]) {
            double theta = (ps[p.F] - ps[y]).angle(ps[x] - ps[y]);
            if (lse(alpha, theta) and lse(theta, beta)) {
                if (ls(cost + p.S, dis[y][p.F]) and ls(cost + p.S, ans)) {
                    dis[y][p.F] = cost + p.S;
                    q.push({-dis[y][p.F], {y, p.F}});
                }
            }
        }
    }
    if (ans > 1e15) ans = -1;
    cout << ans << nl;
}

int32_t main() {
    fast
    const string NAME{"ninja"};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
//        freopen((NAME + ".out").c_str(), "w", stdout);
    }
    setPrec();
    LT solve();
}
