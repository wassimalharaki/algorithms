/**
 * 21:46:33 11/20/24
 * point_3d
 */
// ./algorithms/geometry/OrganizedTemplates/point_3d.cpp
#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;
#define int long long
#define uint unsigned int
#define double long double
#define fast ios_base::sync_with_stdio(false);cin.tie(NULL);cout<<fixed<<setprecision(6);
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

using ftype = double;

//bool eq(const ftype& a, const ftype& b) {
//    return a == b;
//}
bool eq(const ftype& a, const ftype& b) {
    return abs(a - b) < eps;
}

//bool ls(const ftype& a, const ftype& b) {
//    return a < b;
//}
bool ls(const ftype& a, const ftype& b) {
    return (b - a) > eps;
}

bool lse(const ftype& a, const ftype& b) {
    return ls(a, b) or eq(a, b);
}

int sign(const ftype& v) {
    if (ls(0, v)) return 1;
    else if (ls(v, 0)) return -1;
    return 0;
}

ftype sqr(ftype v) { return v * v; }

template<class T>
vec2d<T> matmul(const vec2d<T>& a, const vec2d<T>& b) {
    // Validate dimensions
    if (a.empty() or b.empty() or a[0].size() != b.size()) {
        throw std::invalid_argument("Invalid dimensions for matrix multiplication.");
    }

    size_t rowsA = a.size(), colsA = a[0].size(), colsB = b[0].size();

    vec2d<T> result(rowsA, vec<T>(colsB));

    // Perform matrix multiplication
    for (size_t i = 0; i < rowsA; i++) {
        for (size_t j = 0; j < colsB; j++) {
            for (size_t k = 0; k < colsA; k++) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    return result;
}

struct P {
    ftype x, y, z;
    P(): x{0}, y{0}, z{0} {}
    P(ftype x, ftype y, ftype z): x{x}, y{y}, z{z} {}
    P operator+(const P& p) const {
        return {x + p.x, y + p.y, z + p.z};
    }
    P operator-(const P& p) const {
        return {x - p.x, y - p.y, z - p.z};
    }
    P operator*(ftype v) const {
        return {x * v, y * v, z * v};
    }
    P operator*(const P& p) const { // 3d cross product
        // note that
        // a * b = c
        // where:
        // c.abs() = a.abs() * b.abs() * sin(t)
        // c is orthogonal to the plane formed by a and b
        // where t is the angle between a and b in their common plane.
        P res;
        res.x = y * p.z - z * p.y;
        res.y = z * p.x - x * p.z;
        res.z = x * p.y - y * p.x;
        return res;
    }
    P operator/(ftype v) const {
        return {x / v, y / v, z / v};
    }
    P operator/=(ftype v) {
        return *this = *this / v;
    }
    P& operator*=(ftype v) {
        return *this = *this * v;
    }
    P& operator+=(const P& p) {
        return *this = *this + p;
    }
    P& operator-=(const P& p) {
        return *this = *this - p;
    }
    P& operator*=(const P& p) {
        return *this = *this * p;
    }
    P normalize() const {
        return truncate(1);
    }
    ftype r() const {
        return sqrt(x * x + y * y + z * z);
    }
    ftype dist(const P& p) const {
        return (*this - p).r();
    }
    ftype dot(const P& p) const {
        // note that
        // a.dot(b) = a.abs() * b.abs() * cos(t)
        // where t is the angle between a and b in their common plane.
        return x * p.x + y * p.y + z * p.z;
    }
    ftype xa() const {
        // cos(a) = v.dot(a) / (v.r() * a.r())
        return acos(dot({1, 0, 0})) / r();
    }
    ftype ya() const {
        // cos(a) = v.dot(a) / (v.r() * a.r())
        return acos(dot({0, 1, 0})) / r();
    }
    ftype za() const {
        // cos(a) = v.dot(a) / (v.r() * a.r())
        return acos(dot({0, 0, 1})) / r();
    }
    P applyMat(const vec2d<ftype>& m) const {
        vec2d<ftype> tmp{{x}, {y}, {z}};
        tmp = matmul(m, tmp);
        return {tmp[0][0], tmp[1][0], tmp[2][0]};
    }
    P rotateZccw(ftype a, const P& p = {0, 0, 0}) const {
        P res = *this - p;
        vec2d<ftype> mat{
            {cos(a), -sin(a), 0},
            {sin(a), cos(a), 0},
            {0, 0, 1}
        };

        return res.applyMat(mat) + p;
    }
    P rotateXccw(ftype a, const P& p = {0, 0, 0}) const {
        P res = *this - p;
        vec2d<ftype> mat{
                {1, 0, 0},
                {0, cos(a), -sin(a)},
                {0, sin(a), cos(a)},
        };

        return res.applyMat(mat) + p;
    }
    P rotateYccw(ftype a, const P& p = {0, 0, 0}) const {
        P res = *this - p;
        vec2d<ftype> mat{
                {cos(a), 0, sin(a)},
                {0, 1, 0},
                {-sin(a), 0, cos(a)}
        };

        return res.applyMat(mat) + p;
    }
    P truncate(ftype v) const {
        // returns a vector with norm v and having same direction
        ftype k = r();
        if (sign(k) == 0) return *this;
        v /= k;
        return *this * v;
    }
    bool orthogonal(const P& p) const {
        return eq(dot(p), 0);
    }
};
ostream& operator<<(ostream& out, const P& p) {
    return out << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}
istream& operator>>(istream& in, P& p) {
    int x, y, z;
    in >> x >> y >> z;
    p.x = x;
    p.y = y;
    p.z = z;
    return in;
}

// Ax+By+Cz+D=0
// Where A, B, C are the components of the normal
struct Plane {
    P a, b, c;
    P pvec() const {
        return ((b - a) * (c - a)).normalize();
    }
    Plane(): a{}, b{}, c{} {}
    Plane(const P& a, const P& b, const P& c): a{a}, b{b}, c{c} {}
    bool on(const P& p) const {
        return pvec().orthogonal(p - a);
    }
};
