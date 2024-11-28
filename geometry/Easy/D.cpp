/**
 * 19:00:33 9/4/24
 * D
 */
// ./ICPC/Geometry/Easy/D.cpp
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
#define INF 1000000000000000000ll
#define MOD 1000000007ll
#define PI 3.14159265358979323846264338327950288L
#define EPS 1e-9l
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

void solve() {
    int r, n;
    cin >> r >> n;
    double s = sin(PI / 4);
    cout << r * pow((1 - s) / (1 + s), n) << nl;

    // hyp = sqrt(2 * r^2)
    // C = r + R
    // H = r - R
    // sin(pi / 4) = H / C
    // (r - R) / (r + R) = s
    // r - R = s * r + s * R
    // s * R + R = r - s * r
    // R * (1 + s) = r * (1 - s)
}

int32_t main() {
    fast
    const string NAME{" "};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
        freopen((NAME + ".out").c_str(), "w", stdout);
    }
    LT solve();
}
