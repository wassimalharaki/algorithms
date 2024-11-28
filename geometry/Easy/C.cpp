/**
 * 17:18:33 9/4/24
 * C
 */
// ./ICPC/Geometry/Easy/C.cpp
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

bool ok(int x, const vec<int>& a, int m) {
    int cnt = 0, n = (int)a.size(), i = 0;
    while (i < n - 1 and cnt < m) {
        for (int j = i + 1; j < n; j++) {
            if (a[j] - a[i] >= x) {
                if (cnt == 0) cnt++;
                cnt++;
                i = j;
                break;
            } else if (j == n - 1) {
                i = n - 1;
            }
        }
    }
    return cnt >= m;
}

void solve() {
    int n;
    cin >> n;
    vec<int> a(n);
    for (int& i: a) cin >> i;
    int m;
    cin >> m;
    sort(all(a));
    int z = (1LL << 32) + 1, x = -1;
    for (int b = z; b >= 1; b /= 2) {
        while (ok(x + b, a, m)) x += b;
    }
    cout << x << nl;
}

int32_t main() {
    fast
//    freopen("c.in", "r", stdin);
//    LT
    solve();
}
