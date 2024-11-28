/**
 * 18:03:02 9/10/24
 * E
 */
// ./ICPC/Geometry/Medium/E.cpp
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

void solve() {
    int n, k;
    cin >> n >> k;
    vec<int> p(n);
    for (int& i: p) cin >> i;
    std::sort(p.begin(), p.end());
    int a = 0;
    for (int i = 0; i < n; i++) {
        a = max(a, i + 1 == n ? p[0] - p[i] + 360 : p[i + 1] - p[i]);
    }
    a = 360 - a;
    cout << (a + k - 1) / k << nl;
}

int32_t main() {
    fast
    const string NAME{" "};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
//        freopen((NAME + ".out").c_str(), "w", stdout);
    }
    LT solve();
}
