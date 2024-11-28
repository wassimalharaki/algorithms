#include <bits/stdc++.h>

using namespace std;
#define int long long
#define double long double
#define x real()
#define y imag()

const double eps = 1e-9l, pi = acos(-1.0l);

using ft = double;
using p2 = complex<double>;
using cp = const p2&;
const p2 o = {0 + 0j};

ft dot(cp a, cp b) { return (conj(a) * b).x; }
ft cross(cp a, cp b) { return (conj(a) * b).y; }
ft dist2(cp a, cp b) { return norm(a - b); }
ft diste(cp a, cp b) { return abs(a - b); }
ft angleOfElevation(cp a, cp b) { return arg(b - a); }
ft slope(cp a, cp b) { return tan(angleOfElevation(a, b)); }
p2 toPolar(cp a) { return {abs(a), arg(a)}; }
p2 toCart(cp a) { return polar(a.x, a.y); }
p2 rotateAround(ft t, p2 a, cp r = o) {return (a - r) * polar(1.0l, t) + r;}
// remainder normalizes the angle in [-pi, pi], we get positive non-reflex angle by taking abs.
ft angle(cp a, cp b, cp c) { return abs(remainder(arg(a - b) - arg(c - b), 2.0 * pi)); }
p2 projectVector(cp v, cp p) { return v * dot(p, v) / norm(v); }
p2 projectLine(cp p, cp a, cp b) { return a + (b - a) * dot(p - a, b - a) / norm(b - a); }
p2 reflect(cp p, cp a, cp b) { return a + conj((p - a) / (b - a)) * (b - a); }
p2 intersection(cp a, cp b, cp p, cp q) {
    ft c1 = cross(p - a, b - a), c2 = cross(q - a, b - a);
    return (c1 * q - c2 * p) / (c1 - c2); // undefined if parallel
}
istream& operator>>(istream& in, p2& p) {
    int X, Y;
    in >> X >> Y;
    p = (X, Y);
    return in;
}

void solve() {

}
