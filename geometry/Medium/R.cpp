/**
 * 16:06:39 9/17/24
 * R
 */
// ./ICPC/Geometry/Medium/R.cpp
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
#define P complex<double>
#define X real()
#define Y imag()
#define vec vector
#define LT int T; cin >> T; while (T--)
template<typename T>
using vec2d = vector<vector<T>>;
template<typename T>
using ordered_set = tree<T, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;
using indexed_set = tree<int, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;


double crossProd(P p0, P p1) {
    return (conj(p0) * p1).Y;
}

bool eq(double a, double b) {
    return abs(a - b) <= max(abs(a), abs(b)) * EPS;
}

bool ls(double a, double b) {
    return (b - a) > max(abs(a), abs(b)) * EPS;
}

struct custom_hash {
    static uint64_t splitmix64(uint64_t x) {
        // http://xorshift.di.unimi.it/splitmix64.c
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    }

    size_t operator()(uint64_t x) const {
        static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + FIXED_RANDOM);
    }
};

double floods(const vec<P>& a, int n, unordered_set<int, custom_hash>& flooded) {
    double res = 0, x;
    int i = 0;
    while (i < n) {
        if (not ls(a[i + 1].Y, a[i].Y)) {
            i++;
            continue;
        }
        int st = i;
        double tmp = 0;
        while (i + 1 < n and ls(a[i + 1].Y, a[st].Y)) {
            tmp += crossProd(a[i], a[i + 1]);
            i++;
        }
        i++;
        if (ls(a[st].Y, a[i].Y)) {
            // Compute the intersection point between
            // px(t) = x0 + t * (x1 - x0)
            // py(t) = y0 + t * (y1 - y0)
            // yc = y0 + t * (y1 - y0)
            // (yc - y0) / (y1 - y0) = t
            double t = (a[st].Y - a[i - 1].Y);
            t /= (a[i].Y - a[i - 1].Y);
            x = a[i - 1].X + t * (a[i].X - a[i - 1].X);
            P inter = {x, a[st].Y};
            tmp += crossProd(a[i - 1], inter);
            tmp += crossProd(inter, a[st]);
            res += abs(tmp) / 2;
        } else if (eq(a[i].Y, a[st].Y)) {
            if (flooded.contains(i) or flooded.contains(st)) continue;
            tmp += crossProd(a[i - 1], a[i]);
            tmp += crossProd(a[i], a[st]);
            flooded.insert(n - i - 1);
            flooded.insert(n - st - 1);
            res += abs(tmp) / 2;
        }
    }
    return res;
}

void solve() {
    /*
     * Iterate over the points, when the y of the i + 1 point is less than
     * our y, then we keep on iterating over the points until we reach a y
     * greater than or equal to the y of the original point.
     * calculate the area of the polygon.
     * Continue on from the point with the higher y.
     * Repeat the same process in reverse in relation between: i and i - 1;
     * Make sure to not count the flooded points that were part of the previous pass.
     * This can be done by storing the flooded polygon points
     * (since double counting occurs only in symmetrical cases).
     */
    int n;
    int x, y;
    cin >> n;
    vec<P> a(n);
    for (int i = 0; i < n; i++) {
        cin >> x >> y;
        a[i] = {(double)x, (double)y};
    }
    unordered_set<int, custom_hash> flooded;
    double res = floods(a, n, flooded);
    reverse(a.begin(), a.end());
    res += floods(a, n, flooded);
    cout << res << nl;
}

int32_t main() {
    fast
    const string NAME{" "};
    if (NAME != " ") {
        freopen((NAME + ".in").c_str(), "r", stdin);
    }
    LT solve();
}
