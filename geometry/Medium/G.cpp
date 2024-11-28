/**
 * 21:31:57 9/11/24
 * G
 */
// ./ICPC/Geometry/Medium/G.cpp
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
const int mod = 1e9 + 7;

void solve() {
    const int N = 5010;
    int n;
    cin >> n;
    vec<int> a(N);
    for (int i = 0; i < n; i++) cin >> a[i];
    sort(a.begin(), a.begin() + n);
    vec2d<int> dp(N, vec<int>(N, -1));

    function<int(int, int)> best = [&](int d, int k) -> int {
        if (d == -1) return k == 0;
        int& r = dp[d][k];
        if (r != -1) return r;
        return r = (best(d - 1, k) + best(d - 1, max(0ll, k - a[d]))) % mod;
    };

    int cnt = 0;
    for (int i = 0; i < n; i++) {
        cnt = (cnt + best(i - 1, a[i] + 1)) % mod;
    }
    cout << cnt << nl;
}

int32_t main() {
    fast
    const string NAME{"polygon"};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
    }
    solve();
}
