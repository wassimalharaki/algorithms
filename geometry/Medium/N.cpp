/**
 * 19:11:23 9/7/24
 * N
 */
// ./ICPC/Geometry/Medium/N.cpp
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
const double eps = 1e-9l, pi = 3.14159265358979323846264338327950288l;

bool ok_l(int d, int h, int w, int a) {
    if (d >= 90 * SF) return false;
    double r = d;
    r /= SF;
    r *= pi / 180;
    double ha = a;
    ha /= w;
    ha += w * tan(r) / 2;
    return ha <= h;
}
bool ok_s(int d, int h, int w, int a) {
    double beta = 90 * SF - d;
    beta /= SF;
    beta *= pi / 180;
    double nw = tan(beta) * h;
    double ta = nw * h / 2;
    return ta - a > eps;
}

void solve() {
    int h, w, a;
    cin >> h >> w >> a;
    if (2 * a >= h * w) {
        int x = -1;
        const int z = 90 * SF;
        for (int b = z; b >= 1; b /= 2) {
            while (ok_l(x + b, h, w, a)) x += b;
        }
        double res = x;
        cout << res / SF << nl;
    } else {
        int x = -1;
        const int z = 45 * SF;
        for (int b = z; b >= 1; b /= 2) {
            while (ok_s(x + b, h, w, a)) x += b;
        }
        double res = x;
        cout << res / SF << nl;
    }
}

int32_t main() {
    fast
    const string NAME{"G"};
    freopen((NAME + ".in").c_str(), "r", stdin);
    LT solve();
}
