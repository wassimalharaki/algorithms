/**
 * 13:59:10 11/20/24
 * L
 */
// ./algorithms/geometry/Medium/L.cpp
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
#define allr(v) v.rbegin(), v.rend()
#define NO (cout << "NO" << nl);
#define YES (cout << "YES" << nl);
#define F first
#define S second
#define pb push_back
#define emb emplace_back
#define INF LONG_LONG_MAX
#define MOD 1000000007ll
#define pr pair
#define pii pair<int, int>
#define vec vector
#define LT int T; cin >> T; while (T--)
template<typename T>
using vec2d = vector<vector<T>>;
template<typename T>
using ordered_set = tree<T, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;
using indexed_set = tree<int, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;
const double eps = 1e-9l, pi = acos(-1.0l);

double polyArea(int n, double s) {
    return n * s * s / (4 * tan(pi / n));
}

bool eq(double a, double b) {
    return abs(b - a) <= eps;
}

bool ls(double a, double b) {
    return b - a > eps;
}

bool lse(double a, double b) {
    return ls(a, b) or eq(a, b);
}


void solve() {
    int ip, xr, xd;
    char _;
    cin >> ip >> xr >> _ >> xd;
    double p = ip;
    double x = xr + xd / 1000.0l;
    auto ok = [&](int n) -> bool {
        if (n == 0) return false;
        double a = polyArea(n, p / n);
        return lse(a, p * x);
    };
    int r = -1;
    const int z = 1e6;
    for (int b = z; b >= 1; b >>= 1) {
        while (ok(r + b)) {
            r += b;
            if (r >= z) goto done;
        }
    }
    done:
    if (r <= 0) cout << "Khairy" << nl;
    else if (z <= r) cout << "KEE" << nl;
    else cout << r << nl;
}

int32_t main() {
    fast
    const string NAME{"polygon"};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
//        freopen((NAME + ".out").c_str(), "w", stdout);
    }
    LT solve();
}
