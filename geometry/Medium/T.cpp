/**
 * 16:03:35 9/18/24
 * T
 */
// ./ICPC/Geometry/Medium/T.cpp
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

const double pi = acos(-1.0), eps = 1e-15l;

bool eq(double v0, double v1) {
    return abs(v0 - v1) <= max(v0, v1) * eps;
}

bool ls(double a, double b) {
    return (b - a) > max(abs(a), abs(b)) * eps;
}

double angle(const P& p) {
    double r = atan2(p.Y, p.X);
    if (ls(r, 0)) r += 2 * pi;
    return r;
}

double dif(double a0, double a1) { // a0 - a1
    double a = (a1 - a0);
    a += (ls(pi, a) ? -2 * pi : (ls(a, -pi) ? 2 * pi : 0));
    return a;
}

istream& operator>>(istream& in, P& p) {
    int x, y;
    in >> x >> y;
    p = {(double)x, (double)y};
    return in;
}

int sign(double v) {
    if (ls(v, 0)) return -1;
    else if (ls(0, v)) return 1;
    return 0;
}

void solve() {
    P a, b, c, t;
    int s;
    cin >> a >> b >> c >> t >> s;
    a -= b, c -= b, t -= b;
    double angle_t = angle(t),
           angle_a = angle(a),
           angle_c = angle(c);
    if (not s) { // a can move
        double d = dif(angle_c, angle_t);
        double dac = dif(angle_c, angle_a);
        if (not (ls(0, abs(d)) and ls(abs(d), pi / 2)) or sign(d) != sign(dac)) {
            NO
            return;
        }
        YES
        cout << 180 * (2 * d - dac) / pi << nl;
    } else { // c can move
        double d = dif(angle_a, angle_t);
        double dac = dif(angle_a, angle_c);
        if (not (ls(0, abs(d)) and ls(abs(d), pi / 2)) or sign(d) != sign(dac)) {
            NO
            return;
        }
        YES
        cout << 180 * (2 * d - dac) / pi << nl;
    }
}

int32_t main() {
    fast
    const string NAME{"T"};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
    }
    LT solve();
}
