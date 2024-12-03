#include <bits/stdc++.h>
#define int long long
#define FR(i,a,b) for (int i=(a); i<(b); i++)
#define F(i,n) FR(i,0,n)
using namespace std;
typedef long double ld;
const ld PIl = acos((ld)-1);

typedef int T;
struct p3 {
    T x,y,z;
    p3 operator+(p3 p) {return {x+p.x, y+p.y, z+p.z};}
    p3 operator-(p3 p) {return {x-p.x, y-p.y, z-p.z};}
    p3 operator*(T d) {return {x*d, y*d, z*d};} // mult. by scalar
    T operator|(p3 p) {return x*p.x + y*p.y + z*p.z;} // dot product
    p3 operator*(p3 p) { // cross product
        return {y*p.z - z*p.y, z*p.x - x*p.z, x*p.y - y*p.x};}
};
T sq(p3 p) {return p|p;} // squared norm
ostream& operator<<(ostream& os, p3 p) {
    return os<<"("<<p.x<<","<<p.y<<","<<p.z<<")";}
istream& operator>>(istream& is, p3 &p) {
    return is>>p.x>>p.y>>p.z;}

ld solveBigger(int a, int b, int c, ld disc) {
    assert(disc >= 0);
    if (b > 0)
        return (2*c) / (-b-sqrt(disc));
    else
        return (-b+sqrt(disc)) / (2*a);
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);
    vector<p3> p(2), v(2), w(2);
    F(j,2)
        cin >> p[j] >> v[j] >> w[j];
    p3 d = p[1]-p[0];
    int a = v[0]*v[1]|d, b = (v[0]*w[1] + w[0]*v[1])|d, c = w[0]*w[1]|d;
    if (a < 0) a = -a, b = -b, c = -c;
    int disc = b*b - 4*a*c;
    if (disc < 0) cout<<"never\n";
    else {
        ld cot = solveBigger(a, b, c, disc);
        cout << fixed << setprecision(12);
        cout << PIl / 2 - atan(cot) << "\n";
    }
}