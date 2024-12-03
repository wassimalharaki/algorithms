#include <bits/stdc++.h>
#define FR(i,a,b) for (int i=(a); i<(b); i++)
#define F(i,n) FR(i,0,n)
using namespace std;
typedef long double ld;
const ld PIl = acos((ld)-1);

struct p3 {ld x,y,z;};
p3 operator+(p3 a, p3 b) {return {a.x+b.x, a.y+b.y, a.z+b.z};}
p3 operator-(p3 a, p3 b) {return {a.x-b.x, a.y-b.y, a.z-b.z};}
p3 operator*(p3 a, ld d) {return {a.x*d, a.y*d, a.z*d};}

ld eval(function<ld(ld)> f, ld a, ld b) {
    return (b-a)/6 * (f(a) + 4*f((a+b)/2) + f(b));
}

ld integrate(function<ld(ld)> f, ld a, ld b) {
    ld mid = (a+b)/2, q1 = eval(f,a,b),
            q2 = eval(f,a,mid) + eval(f,mid,b);
    if (abs(q2-q1) <= 1e-11)
        return q2 + (q2-q1)/15;
    return integrate(f,a,mid) + integrate(f,mid,b);
}

ld cover(p3 a, p3 b) {
    ld cross = a.x*b.y - a.y*b.x;
    auto f = [&](ld t) {
        p3 p = a+(b-a)*t;
        ld sqXY = p.x*p.x + p.y*p.y;
        ld r = .5 - atan(p.z/sqrt(sqXY))/PIl;
        ld dphi = cross/sqXY;
        return r*r*dphi/2;
    };
    ld res = integrate(f, 0, 1);
    return res;
}

p3 sph(ld lat, ld lon) {
    lat *= PIl/180, lon *= PIl/180;
    return {cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)};
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);
    int n; cin >> n;
    vector<p3> ps(n);
    F(i,n) {
        int lat, lon; cin >> lat >> lon;
        ps[i] = sph(lat, lon);
    }
    ld sum = 0;
    F(i,n) sum += cover(ps[i], ps[(i+1)%n]);
    cout << fixed << setprecision(12) << sum << "\n";
}