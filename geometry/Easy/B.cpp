/**
 * 16:59:11 9/4/24
 * B
 */
// ./ICPC/Geometry/Easy/B.cpp
#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;
#define int long long
#define uint unsigned int
#define double long double
#define fast ios_base::sync_with_stdio(false);cin.tie(NULL);
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
#define P pair<double, double>
#define X F
#define Y S
#define vec vector
#define LT int T; cin >> T; while (T--)
template<typename T>
using vec2d = vector<vector<T>>;
template<typename T>
using ordered_set = tree<T, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;
using indexed_set = tree<int, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;

P pyth(double a, double hyp) {
    a /= 180;
    a *= 3.14159265358979323846264338327950288L;
    double s = sin(a), c = cos(a);
    return {hyp * c, hyp * s};
}

void solve() {
    int n;
    cin >> n;
    vec<P> a(n);
    for (auto& p: a) cin >> p.X >> p.Y;
    P s{0, 0};
    for (auto& p: a) {
        auto d = pyth(p.F, p.S);
        s.X += d.F;
        s.Y += d.S;
    }
    cout << setprecision(32) << sqrt(s.X * s.X + s.Y * s.Y) << nl;
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
