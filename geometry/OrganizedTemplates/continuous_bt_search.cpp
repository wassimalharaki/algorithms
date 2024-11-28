/**
 * 21:31:51 9/25/24
 * binarySearch
 */
// ./ICPC/Geometry/Templates/binarySearch.cpp
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

/*
 * Ternary search

double r_lo = 0;    ll lo = reinterpret_cast<ll&>(r_lo);
double r_hi = MAXV; ll hi = reinterpret_cast<ll&>(r_hi);
while (hi > lo)
{
    ll lm = lo + (hi - lo)/3;
    ll rm = hi - (hi - lo)/3;
    if (calc(reinterpret_cast<double&>(lm)) < calc(reinterpret_cast<double&>(rm))) hi = rm-1;
    else lo = lm+1;
}
 */

using ui = uint64_t;
#define cui(v) reinterpret_cast<ui&>(v)
#define cd(v) reinterpret_cast<double&>(v)

void solveP() {
    // Note that in case of negative numbers, we need to convert from two's complement to keep order.
    double l_real = 0, r_real = 1e15l;
    ui l = cui(l_real), r = cui(r_real);
    auto check = [&](double v) -> bool {
        return v == 0; // Add logic here
    };
    while (l < r) {
        ui mid = l + (r - l) / 2;
        if (check(cd(mid))) {
            r = mid;
        } else {
            l = mid + 1;
        }
    }
}

const double eps = 1e-9l;
bool ls(double a, double b) {
    return (b - a) > max(abs(a), abs(b)) * eps;
}

//double ternary_search(double l, double r) {
//    while (r - l > eps) {
//        double m1 = l + (r - l) / 3;
//        double m2 = r - (r - l) / 3;
//        double f1 = f(m1);      //evaluates the function at m1
//        double f2 = f(m2);      //evaluates the function at m2
//        if (f1 < f2)
//            l = m1;
//        else
//            r = m2;
//    }
//    return f(l);                    //return the maximum of f(x) in [l, r]
//}


// Where LOG is the iterations
// const int LOG = 200;
//for(int i = 0; i < LOG; i++)
//{
//long double m1 = (A * 2 + B) / 3.0;
//long double m2 = (A + 2 * B) / 3.0;
//
//if(f(m1) > f(m2))
//A = m1;
//else
//B = m2;
//}
//
//ans = f(A);

//while(B - A > 4)
//{
//int m1 = (A + B) / 2;
//int m2 = (A + B) / 2 + 1;
//
//if(f(m1) > f(m2))
//A = m1;
//else
//B = m2;
//}
//
//ans = inf;
//
//for(int i = A; i <= B; i++)
//ans = min(ans , f(i));

void solve() {
    double c;
    cin >> c;
    double l = 0, r = 1e15l;
    const int iter = 400; // specify
    auto ok = [&](double v) -> bool {
        return not ls(v * v + sqrt(v) - c, 0);
    };
    for (int i = 0; i < iter and ls(l, r); i++) {
        double m = (l + r) * 0.5l;
        if (ok(m)) r = m;
        else l = m;
    }
    // or (l + r) * 0.5;
    cout << (l + r) * 0.5l << nl;
}
