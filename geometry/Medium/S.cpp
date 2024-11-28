/**
 * 16:08:05 9/13/24
 * S
 */
// ./ICPC/Geometry/Medium/S.cpp
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
#define P complex<int>
#define X real()
#define Y imag()
#define vec vector
#define LT int T; cin >> T; while (T--)
template<typename T>
using vec2d = vector<vector<T>>;
template<typename T>
using ordered_set = tree<T, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;
using indexed_set = tree<int, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;

const double pi = 3.14159265358979323846264338327950288l, eps = 1e-6;

double sqr(double v) { return v * v; }
bool eq(double v0, double v1) {
    return abs(v0 - v1) <= eps;
}


void solve() {
    int x0, y0, r0, x1, y1, r1;
    cin >> x0 >> y0 >> r0 >> x1 >> y1 >> r1;
    double d = sqrt(sqr(x0 - x1) + sqr(y0 - y1));
    if ((x0 == x1 and y0 == y1) or d >= r0 + r1) {
        fail:
        cout << -1 << nl;
        return;
    } else if (d <= abs(r0 - r1)) {
        if (not eq(min(r1, r0) + d, max(r0, r1))) goto fail;
    }
    double a = sqr(r0) * acos((sqr(d) + sqr(r0) - sqr(r1)) / (2 * d * r0))
            + sqr(r1) * acos((sqr(d) + sqr(r1) - sqr(r0)) / (2 * d * r1))
            - sqrt((-d + r0 + r1) * (d + r0 - r1) * (d - r0 + r1) * (d + r0 + r1)) / 2;

    cout << sqr(max(r0, r1)) * pi - a << nl;
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
