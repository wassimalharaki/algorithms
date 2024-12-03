#include <bits/stdc++.h>
#define FR(i,a,b) for (int i=(a); i<(b); i++)
#define F(i,n) FR(i,0,n)
using namespace std;
const int INF = 1e9, T=3;
typedef long double ld;

struct p3 {ld x,y,z;};
p3 operator+(p3 a, p3 b) {return {a.x+b.x, a.y+b.y, a.z+b.z};}
p3 operator-(p3 a, p3 b) {return {a.x-b.x, a.y-b.y, a.z-b.z};}
p3 operator*(p3 a, ld d) {return {a.x*d, a.y*d, a.z*d};}
p3 operator/(p3 a, ld d) {return a*(1/d);}
bool operator<(p3 a, p3 b) {return tie(a.x, a.y, a.z) < tie(b.x, b.y, b.z);}
istream& operator>>(istream& is, p3 &p) {
    return is>>p.x>>p.y>>p.z;}
ostream& operator<<(ostream& os, p3 p) {
    return os<<"("<<p.x<<","<<p.y<<","<<p.z<<")";}
p3 operator*(p3 a, p3 b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
ld operator|(p3 a, p3 b) {return a.x*b.x + a.y*b.y + a.z*b.z;}
ld abs(p3 a) {return sqrt(a|a);}

struct edge {int v,lo,hi;};
struct uh {int u,h;};
bool operator<(uh a, uh b) {return a.h < b.h;};

p3 planeInter(p3 a, p3 b, int h) {
    return (a*(b.z-h) - b*(a.z-h))/(b.z-a.z);
}

vector<p3> polygonUnder(vector<p3> ps, int h) {
    vector<p3> under;
    for (int i=0, n=ps.size(); i<n; i++) {
        p3 a = ps[i], b = ps[(i+1)%n];

        // Keep points under h
        if (a.z <= h)
            under.push_back(a);

        // Add intersections
        if ((a.z-h)*(b.z-h) < 0) // if a and b on either side
            under.push_back(planeInter(a,b,h));
    }
    return under;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);
    int n; cin>>n;
    map<pair<p3,p3>, int> fromEdge;
    vector<vector<p3>> ps(n, vector<p3>(T));
    vector<vector<edge>> g(n);
    F(i,n) {
        F(j,T) cin>>ps[i][j];
        F(j,T) {
            p3 a = ps[i][j], b = ps[i][(j+1)%T];
            if (fromEdge.count({a,b})) {
                int k = fromEdge[{a,b}];
                int lo = min(a.z, b.z), hi = max(a.z, b.z);
                g[i].push_back({k,lo,hi});
                g[k].push_back({i,lo,hi});
            } else {
                fromEdge[{b,a}] = i;
            }
        }
    }
    vector<int> h(n, -INF);
    h[0] = ps[0][0].z;
    priority_queue<uh> pq;
    pq.push({0,h[0]});
    while (!pq.empty()) {
        uh cur = pq.top(); pq.pop();
        int u = cur.u;
        if (cur.h < h[u])
            continue;
        for (edge e : g[u]) if (h[u] >= e.lo) {
                int newH = min(h[u], e.hi);
                if (h[e.v] < newH) {
                    h[e.v] = newH;
                    pq.push({e.v, newH});
                }
            }
    }
    ld tot = 0;
    F(i,n) {
        vector<p3> under = polygonUnder(ps[i], h[i]);
        p3 sum{0,0,0};
        int m = under.size();
        F(j,m)
            sum = sum + under[j]*under[(j+1)%m];
        tot += abs(sum)/2;
    }
    cout<<fixed<<setprecision(12)<<tot<<"\n";
}