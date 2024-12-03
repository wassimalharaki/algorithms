#include <bits/stdc++.h>
#define FR(i,a,b) for (int i=(a); i<(b); i++)
#define F(i,n) FR(i,0,n)
using namespace std;
typedef long double ld;

typedef complex<ld> pt;
ld cross(pt a, pt b) {return imag(conj(a)*b);}

// area of [x,oo) x [y,oo) inter C((0,0),r)
ld quadCircle(ld r, ld x, ld y) {
    assert(x >= 0 && y >= 0);
    if (x*x+y*y >= r*r)
        return 0;
    pt p(x,y), a(sqrt(r*r-y*y), y), b(x, sqrt(r*r-x*x));
    return (cross(p,a) + r*r*asin(min(cross(a,b)/r/r, 1.0L)) + cross(b,p))/2;
}

// area of [x,x+1] x [y,y+1] inter C((0,0),r)
ld squareCircle(ld r, ld x, ld y) {
    return quadCircle(r,x,y) - quadCircle(r,x+1,y) - quadCircle(r,x,y+1) + quadCircle(r,x+1,y+1);
}

struct p3 {int x,y,z;};

ld getRadius(ld d, ld ratio) {
    ld alpha = acos(min(cbrt(d*d*ratio), 1.0L));
    return d*tan(alpha);
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);
    int n,l; cin>>n>>l;
    vector<p3> ps(n);
    for (p3 &p : ps)
        cin >> p.x >> p.y >> p.z;
    ld lo = 0, hi = 1; // ratio 3D angle / paint used
    for (;;) {
        ld mid = (lo+hi)/2;
        vector<ld> paintOn(n,0);
        F(i,n) {
            p3 p = ps[i];
            auto paintFace = [&](int d, int x, int y) {
                if (d > 0)
                    paintOn[i] += squareCircle(getRadius(d,mid), x, y);
            };
            paintFace(p.x, p.y, p.z);
            paintFace(p.y, p.x, p.z);
            paintFace(p.z, p.x, p.y);
        }
        if (mid == hi || mid == lo) {
            cout << fixed << setprecision(8);
            F(i,n) cout << paintOn[i] << "\n";
            break;
        } else {
            ld sum = 0;
            F(i,n) sum += paintOn[i];
            if (sum < l)
                hi = mid;
            else
                lo = mid;
        }
    }
}