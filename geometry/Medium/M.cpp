/**
 * 16:05:22 11/20/24
 * M
 */
// ./algorithms/geometry/Medium/M.cpp
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

struct DSU {
    vector<int> p;

    DSU(int n) { p.resize(n, -1); }

    int find(int x) {
        return p[x] < 0 ? x : p[x] = find(p[x]);
    }

    int size(int x) { return - p[find(x)]; }

    bool same_set(int x, int y) {
        return find(x) == find(y);
    }

    bool merge(int x, int y) {
        x = find(x); y = find(y);
        if (x == y) return false;
        if (p[x] > p[y]) swap(x, y);
        p[x] += p[y]; p[y] = x;
        return true;
    }
};

struct P {
    int x, y, z;
    P(int x, int y, int z): x{x}, y{y}, z{z} {}
    P(): P(0, 0, 0) {}
    bool operator<(const P& p) const {
        if (x != p.x) return x < p.x;
        if (y != p.y) return y < p.y;
        return z < p.z;
    }
    bool operator>(const P& p) const {
        if (x != p.x) return x > p.x;
        if (y != p.y) return y > p.y;
        return z > p.z;
    }
    bool operator==(const P& p) const {
        return x == p.x and y == p.y and z == p.z;
    }
};
istream& operator>>(istream& in, P& p) {
    return in >> p.x >> p.y >> p.z;
}
struct S {
    P a, b;
    S(P a, P b) {
        if (a > b) swap(a, b);
       this->a = a;
       this->b = b;
    }
    bool operator<(const S& s) const {
        if (a != s.a) return a < s.a;
        return b < s.b;
    }
};


void solve() {
    int n;
    cin >> n;
    map<S, int> m;
    DSU dsu(n);
    auto add = [&](const S& s, int i) -> void {
        if (m.find(s) == m.end()) m[s] = i;
        else dsu.merge(m[s], i);
    };

    for (int i = 0; i < n; i++) {
        P a, b, c;
        cin >> a >> b >> c;
        S ab(a, b), bc(b, c), ca(c, a);
        add(ab, i);
        add(bc, i);
        add(ca, i);
    }
    int sz = dsu.size(0);
    if (sz == n) YES
    else NO
}

int32_t main() {
    fast
    const string NAME{"triangles"};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
//        freopen((NAME + ".out").c_str(), "w", stdout);
    }
    LT solve();
}
