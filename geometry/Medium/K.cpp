/**
 * 21:35:33 9/11/24
 * K
 */
// ./ICPC/Geometry/Medium/K.cpp
#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;
#define int long long
#define uint unsigned int
#define double long double
#define fast ios_base::sync_with_stdio(false);cin.tie(NULL);cout<<setprecision(25);
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
#define P complex<double>
#define X real()
#define Y imag()
#define vec vector
#define LT int T; cin >> T; while (T--)
template<typename T>
using vec2d = vector<vector<T>>;
template<typename T>
using ordered_set = tree<T, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;
using indexed_set = tree<int, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;

const double pi = acos(-1.0);
const double eps = 1e-7l;
P rotate(const P& p0, double a) {
    double xs = p0.X, ys = p0.Y;
    return {xs * cos(a) - ys * sin(a), xs * sin(a) + ys * cos(a)};
}

double angle(const P& p) {
    double r = atan2(p.Y, p.X);
    if (r < 0) r += 2 * pi;
    return r;
}

double sqr(double v) { return v * v; }

double dot(const P& p0, const P& p1) {
    return p0.X * p1.X + p0.Y * p1.Y;
}

bool eq(double v0, double v1) {
    return abs(v0 - v1) <= max(v0, v1) * eps;
}

bool leq(double a, double b) {
    return (b - a) >= max(abs(a), abs(b)) * eps;
}

bool eq(const P& p0, const P& p1) {
    return eq(p0.X, p1.X) and eq(p0.Y, p1.Y);
}

vec<P> intersection(const pair<P, P>& s, double r, P c) {
    P e = s.F, l = s.S;
    P d = l - e, f = e - c;
    auto gp = [&](double t) -> P { return e + t * d; };
    double A = dot(d, d),
           B = 2 * dot(f, d),
           C = dot(f, f) - r * r;

    double delta = B * B - 4 * A * C;
    if (delta < 0) return {};
    else if (eq(delta, 0)) return {gp(-B / (2 * A))};
    delta = sqrt(delta);
    double t0 = (-B - delta) / (2 * A),
           t1 = (-B + delta) / (2 * A);
    vec<P> res;
    if (leq(0, t0) and leq(t0, 1)) res.emplace_back(gp(t0));
    if (leq(0, t1) and leq(t1, 1)) res.emplace_back(gp(t1));
    return res;
}

double len(const pair<P, P>& s) {
    const auto& p0 = s.F,
                p1 = s.S;
    double x = p0.X - p1.X, y = p0.Y - p1.Y;
    return sqrt(x * x + y * y);
}

void counterclock(double& an, double a) {
    if (leq(2 * pi, an)) an -= 2 * pi;
    if (leq(a, an)) a += 2 * pi;
    an = a - an;
}

void clock(double& an, double a) {
    if (leq(an, 0)) an += 2 * pi;
    if (leq(a, an)) an = an - a;
    else an = 2 * pi - (a - an);
}

double calc(bool isa, double an, double r, double a, double x, double y, double ia, double ib) {
    P p = rotate({r, 0}, a), sp = {x, y};
    pair<P, P> v = {sp, p};
    auto in = intersection(v, r, {0, 0});
    vec<P> tmp;
    for (P& ptmp: in) {
        if (not eq(ptmp, p)) {
            tmp.push_back(ptmp);
        }
    }
    if (tmp.empty()) { // only intersection point is p itself
        return len(v) + r;
    } else {
        P p0 = tmp[0];
        double anp0 = angle(p0);
        if (leq(ia, anp0) and leq(anp0, ib)) {
            return len(v) + r;
        }
    }
    double d = hypot(x, y);
    double tand = sqrt(sqr(d) - sqr(r));
    double rotan = acos(r / d);
    double oan = an;
    if (isa) {
        an += rotan;
        counterclock(an, a);
        counterclock(oan, a);
    } else {
        an -= rotan;
        clock(an, a);
        clock(oan, a);
    }
    return min(an * r + tand + r, oan * r + d);
}

void solve() {
    int r, x, y, a, b;
    cin >> r >> a >> b >> x >> y;
    double rf, xf, yf, af, bf;
    // Store as floats
    rf = r, xf = x, yf = y, af = a, bf = b;

    // Change to radians
    af *= pi / 180, bf *= pi / 180;

    // Calculate the counter-closkwise angle formed by (x, y) with x-axis.
    double an = angle({xf, yf});

    // Calculate the distance from the point to the center
    double d = hypot(xf, yf);

    // If the point's angle is between alpha and beta, then go straight to the center
    if ((leq(af, an) and leq(an, bf)) or d < r) {
        cout << d << nl;
        return;
    }

    // Calculate the distance from the point to the starting points of the gap and take the minimum of which
    double bd = calc(false, an, rf, bf, xf, yf, af, bf);
    double ad = calc(true, an, rf, af, xf, yf, af, bf);
    assert(bd > 0);
    assert(ad > 0);
    assert(bd == bd);
    assert(ad == ad);
    cout << (leq(ad, bd) ? ad : bd) << nl;
}


int32_t main() {
    fast
    const string NAME{"walls"};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
//        freopen((NAME + ".out").c_str(), "w", stdout);
    }
    LT solve();
}

