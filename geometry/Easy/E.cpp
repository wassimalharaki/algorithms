/**
 * 19:25:56 9/4/24
 * E
 */
// ./ICPC/Geometry/Easy/E.cpp
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

int c2(int n) {
    if (n < 2) return 0;
    return (n * (n - 1)) >> 1;
}

void solve() {
    int n;
    cin >> n;
    vec<int> a(n);
    for (int& i: a) cin >> i;
    map<int, int> m;
    for (int i: a) m[i]++;
    vec<int> c;
    for (auto& p: m) c.push_back(c2(p.S));
    int sz = (int)c.size();
    vec<int> pre(sz + 1);
    for (int i = 0; i < sz; i++) pre[i + 1] = pre[i] + c[i];
    int r = 0;
    for (int i = 0; i < sz - 1; i++) {
        r += c[i] * (pre[sz] - pre[i + 1]);
    }
    cout << r << nl;
}

int32_t main() {
    fast
    freopen("brimore.in", "r", stdin);
    solve();
}
