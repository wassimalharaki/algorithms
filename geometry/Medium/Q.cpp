/**
 * 21:22:28 9/10/24
 * Q
 */
// ./ICPC/Geometry/Medium/Q.cpp
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

const int SF = 1e9;
const double eps = 1e-10l;

bool ok(int pos, int h, int v, int t) {
    if (pos > h * SF) return false;
    double p = pos;
    p /= SF;
    double hyp = hypot(p, v);
    return hyp <= t;
}

void solve() {
    int h, v, t;
    cin >> h >> v >> t;
    const int z = h * SF;
    int x = -1;
    for (int b = z; b >= 1; b /= 2) {
        while (ok(x + b, h, v, t)) x += b;
    }
    double a = x;
    a /= SF;
    double d = hypot(a, (double)v) + h - a;
    cout << d << nl;
}

int32_t main() {
    fast
    freopen("cross.in", "r", stdin);
    LT solve();
}
