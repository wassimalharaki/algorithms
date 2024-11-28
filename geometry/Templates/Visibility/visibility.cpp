///**
// /tmp"
// * 10:43:50 11/11/24
// * point_vis
// */
//// ./ICPC/Geometry/OrganizedTemplates/point_vis.cpp
////#include "Point.cpp"
///**
// * 15:15:28 9/27/24
// * Point
// */
//// ./ICPC/Geometry/OrganizedTemplates/Point.cpp
//#include <bits/stdc++.h>
//#include <ext/pb_ds/assoc_container.hpp>
//#include <ext/pb_ds/tree_policy.hpp>
//
//using namespace std;
//using namespace __gnu_pbds;
//#define int long long
//#define uint unsigned int
//#define double long double
//#define fast ios_base::sync_with_stdio(false);cin.tie(NULL);
//#define nl '\n'
//#define all(v) v.begin(), v.end()
//#define NO (cout << "NO" << nl);
//#define YES (cout << "YES" << nl);
//#define F first
//#define S second
//#define INF LONG_LONG_MAX
//#define MOD 1000000007ll
//#define EPS 1e-9l
//#define PI 3.14159265358979323846264338327950288L
//#define pii pair<int, int>
//#define X real()
//#define Y imag()
//#define vec vector
//#define LT int T; cin >> T; while (T--)
//template<typename T>
//using vec2d = vector<vector<T>>;
//template<typename T>
//using ordered_set = tree<T, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;
//using indexed_set = tree<int, null_type, less<>, rb_tree_tag, tree_order_statistics_node_update>;
//
//namespace Geometry {
//    using ftype = double;
//    using atype = double;
//    const double eps = 1e-9l;
//    const double pi = acos(-1.0l);
//    const double e = exp(1.0L);
//    const double inf = INF;
//    const int iinf = LONG_LONG_MAX;
//    const int prec = 25;
//
//    ftype sqr(ftype v) { return v * v; }
//
//    void setPrec(ostream& out = cout) {
//        out << fixed << setprecision(prec);
//    }
//
//    //bool eq(const ftype& a, const ftype& b) {
//    //    return a == b;
//    //}
//    bool eq(const ftype& a, const ftype& b) {
//        return abs(a - b) < eps;
//    }
//
//    //bool ls(const ftype& a, const ftype& b) {
//    //    return a < b;
//    //}
//    bool ls(const ftype& a, const ftype& b) {
//        return (b - a) > eps;
//    }
//
//    //bool lse(const ftype& a, const ftype& b) {
//    //    return a <= b;
//    //}
//    bool lse(const ftype& a, const ftype& b) {
//        return ls(a, b) or eq(a, b);
//    }
//
//    //bool gt(const ftype& a, const ftype& b) {
//    //    return a > b;
//    //}
//    bool gt(const ftype& a, const ftype& b) {
//        return not ls(a, b) and not eq(a, b);
//    }
//
//    // Returns angle with respect to the x-axis.
//    atype angle(atype x, atype y) {
//        return atan2(y, x);
//    }
//    atype torad(atype a) { return a / 180 * pi; }
//    atype todeg(atype a) { return a / pi * 180; }
//
//    int sign(const ftype& v) {
//        if (ls(0, v)) return 1;
//        else if (ls(v, 0)) return -1;
//        return 0;
//    }
//
//    // Orientation of 3 points forming an angle
//    enum class Orientation {
//        l = 1,
//        r = -1,
//        c = 0
//    };
//
//    // 2D point
//    struct P {
//        ftype x, y;
//        static P polar(const ftype& r, const atype& a) {
//            return {r * cos(a), r * sin(a)};
//        }
//        P(): x{0}, y{0} {}
//        P(ftype x, ftype y): x{x}, y{y} {}
//        explicit P(const complex<ftype>& p): P(p.real(), p.imag()) {}
//        explicit P(const pair<ftype, ftype>& p): P(p.F, p.S) {}
//        P(const P& p): P(p.x, p.y) {}
//        explicit operator complex<ftype>() const {
//            return {x, y};
//        }
//        explicit operator pair<ftype, ftype>() const {
//            return {x, y};
//        }
//        ftype dist(const P& p) const {
//            return (*this - p).abs();
//        }
//        atype r() const {
//            return hypot(x, y);
//        }
//        atype a() const {
//            return angle(x, y);
//        }
//        atype ua(const P& p) const {
//            // undirected angle
//            double v0 = fabs(a() - p.a()), v1 = fabs(p.ap() - ap());
//            if (ls(v0, v1)) return v0;
//            else return v1;
//        }
//        atype ap() const {
//            atype res = angle(x, y);
//            if (res < 0) res += 2 * pi;
//            return res;
//        }
//        P operator-() const {
//            return *this * -1;
//        }
//        P operator+() const {
//            return *this;
//        }
//        P operator+=(const P& p) {
//            x += p.x;
//            y += p.y;
//            return *this;
//        }
//        P operator-=(const P& p) {
//            x -= p.x;
//            y -= p.y;
//            return *this;
//        }
//        P operator*=(const ftype& v) {
//            x *= v;
//            y *= v;
//            return *this;
//        }
//        P operator/=(const ftype& v) {
//            x /= v;
//            y /= v;
//            return *this;
//        }
//        P operator%=(const ftype& v) {
//            x = fmod(x, v);
//            y = fmod(y, v);
//            return *this;
//        }
//        P operator^=(const ftype& an) {
//            return *this = rotateccw(an);
//        }
//        P operator|=(const ftype& an) {
//            return *this = rotatecw(an);
//        }
//        P operator|(const ftype& an) {
//            P res = *this;
//            res |= an;
//            return res;
//        }
//        P operator^(const ftype& an) {
//            P res = *this;
//            res ^= an;
//            return res;
//        }
//        P operator+(const P& p) const {
//            P res = *this;
//            res += p;
//            return res;
//        }
//        P operator-(const P& p) const {
//            P res = *this;
//            res -= p;
//            return res;
//        }
//        P operator*(const ftype& v) const {
//            P res = *this;
//            res *= v;
//            return res;
//        }
//        P operator/(const ftype& v) const {
//            P res = *this;
//            res /= v;
//            return res;
//        }
//        P operator%(const ftype& v) const {
//            P res = *this;
//            res %= v;
//            return res;
//        }
//        bool operator==(const P& p) const {
//            return eq(x, p.x) and eq(y, p.y);
//        }
//        bool operator<(const P& p) const {
//            return eq(p.x, x) ? ls(y, p.y) : ls(x, p.x);
//        }
//        bool operator<=(const P& p) const {
//            return *this < p or *this == p;
//        }
//        ftype cross(const P& b, const P& c) const {
//            return (b - *this).cross(c - *this);
//        }
//        Orientation orientation(const P& pb, const P& pc) const {
//            const P& pa = *this;
//            ftype d = (pb - pa).cross(pc - pa);
//            return static_cast<Orientation>(
//                    ls(0, d) - ls(d, 0)
//            );
//        }
//        ftype cross(const P& p) const {
//            return x * p.y - y * p.x;
//        }
//        ftype dot(const P& p) const {
//            return x * p.x + y * p.y;
//        }
//        ftype norm() const {
//            return dot(*this);
//        }
//        atype abs() const {
//            return sqrt((atype)norm());
//        }
//        P truncate(ftype v) const {
//            // returns a vector with norm v and having same direction
//            ftype k = abs();
//            if (sign(k) == 0) return *this;
//            v /= k;
//            return {x * v, y * v};
//        }
//        P normalize() const {
//            return truncate(1);
//        }
//        P normal_l() const {
//            return {-y, x};
//        }
//        P normal_r() const {
//            return {y, -x};
//        }
//        P rotateccw(const atype& angle, const P& ref = {0, 0}) const {
//            P res;
//            ftype xs = x - ref.x, ys = y - ref.y;
//            res.x = ref.x + (xs * cos(angle) - ys * sin(angle));
//            res.y = ref.y + (xs * sin(angle) + ys * cos(angle));
//            return res;
//        }
//        P rotatecw(const atype& angle, const P& ref = {0, 0}) const {
//            return rotateccw(-angle, ref);
//        }
//        friend P operator*(const ftype& v, const P& p) {
//            return p * v;
//        }
//        friend istream& operator>>(istream& in, P& p) {
//            int xv, yv;
//            in >> xv >> yv;
//            p = P(xv, yv);
//            return in;
//        }
//        friend ostream& operator<<(ostream& out, const P& p) {
//            setPrec(out);
//            return out << "(" << p.x << ", " << p.y << ")";
//        }
//    };
//
//    // basic comparison functions
//    auto angleCmp = [](const P& p0, const P& p1) -> bool {
//        return ls(p0.a(), p1.a());
//    };
//    auto radiusCmp = [](const P& p0, const P& p1) -> bool {
//        return ls(p0.r(), p1.r());
//    };
//    auto stdCmp = [](const P& p0, const P& p1) -> bool {
//        return p0 < p1;
//    };
//    auto xCmp = [](const P& p0, const P& p1) -> bool {
//        return ls(p0.x, p1.x);
//    };
//    auto yCmp = [](const P& p0, const P& p1) -> bool {
//        return ls(p0.y, p1.y);
//    };
//
//    // Hash function
//    struct p_hash {
//        static __int128 splitmix128(__int128 x) {
//            // gr = 0x9e3779b97f4a7c15f39cc0605cedc835
//            // c0 = 0xbf58476d1ce4e5b9a3f7b72c1e3c9e3b
//            // c1 = 0x94d049bb133111ebb4093822299f31d7
//            __int128 gr = 0x9e3779b97f4a7c15, c0 = 0xbf58476d1ce4e5b9, c1 = 0x94d049bb133111eb;
//            gr <<= 64;
//            c0 <<= 64;
//            c1 <<= 64;
//            gr += 0xf39cc0605cedc835;
//            c0 += 0xa3f7b72c1e3c9e3b;
//            c1 += 0xb4093822299f31d7;
//            x += gr;
//            x = (x ^ (x >> 62)) * c0;
//            x = (x ^ (x >> 59)) * c1;
//            return x ^ (x >> 63);
//        }
//
//        size_t operator()(const P& p) const {
//            static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
//            __int128 x = (((__int128)p.x) << 64) | ((__int128)p.y);
//            return splitmix128(x + FIXED_RANDOM);
//        }
//    };
//}
//using namespace Geometry;
//
////如果两个向量的叉积结果符号相反（一个为正，一个为负），则说明线段AC和线段AD在线段AB的两侧，即线段AC和线段AD相交。
//double RayIntersect(const P& a, const P& b, const P& c, const P& d, int* sides = NULL) {
//    double cp1 = (c-a).cross(b-a), cp2 = (d-a).cross(b-a);
//    double dp1 = (c-a).dot(b-a), dp2 = (d-a).dot(b-a);
//    if (sides) *sides = (cp1 < -eps || cp2 < -eps) + 2 * (cp1 > eps || cp2 > eps);//sides 0,1,2,3,如果两个线段相交,一定是3  1,2是ab直线穿过c或者d的情况
//    if (cp1 < -eps && cp2 < -eps || cp1 > eps && cp2 > eps) return -1.0;//叉积结果同号，线段不相交
//    return (abs(cp1) < eps && abs(cp2) < eps) ? 0 : (dp1*cp2-dp2*cp1)/(cp2-cp1);//a到交点的距离乘以ab长度
//}
//
//bool POnLine(const P& a, const P& b, const P& p) {
//    double ln = (b-a).abs(), cp = (b-a).cross(p-a), dp = (b-a).dot(p-a);
//    return abs(cp/ln) < eps && dp/ln > -eps && dp/ln < ln+eps;//cp=0 && 0<dp<ln^2
//}
//
//int32_t main(){
//    freopen("A.in", "r", stdin);
//    int n;
//    cin>>n;
//    vector<P> polygon(n);
//    P s,g;
//    for (int i = 0; i < n; ++i) {
//        cin>>polygon[i];
//    }
//    cin>>g>>s;
//    vector<P> endpoints;
//    //找指向每个顶点方向距离最近的intersecting边
//    polygon.push_back(g);//g点方向也需要找
//    for(auto p:polygon){
//        vector<tuple<double,int>> intersections;//交点距离，遮挡面
//        if((p-s).abs()<eps)continue;
//        p=(p-s)/(p-s).abs()+s;
//        //找出sculpture向各个顶点出发, 不被遮挡能走的最远距离maxd，即视野不被遮挡的区域的终点(在多边形的边上)
//        for(int i=0;i<n;++i){
//            int side=0;
//            auto d= RayIntersect(s,p,polygon[i],polygon[(i+1)%n],&side);
//            if(d>-eps)
//                intersections.emplace_back(d,side);
//        }
//        ranges::sort(intersections);
//        double maxd=0;
//        int sides=0;
//        for(auto [d,s]:intersections){
//            maxd=d;
//            sides|=s;
//            if(sides==3)break;
//        }
//        endpoints.push_back(s+(p-s)*maxd);//从sculpture出发,沿着射线,视野不被遮挡的区域的终点(在多边形的边上)
//        //计算顶点p2，到sculpture发出的穿过顶点i的射线，的距离,以及交点
//        for(auto p2:polygon){
//            auto ortho=P(s.y-p.y,p.x-s.x)*1e5;
//            auto d= RayIntersect(s,p,p2,p2+ortho);
//            if(d<-eps)
//                d= RayIntersect(s,p,p2,p2-ortho);
//            if(d>eps && d<maxd-eps)
//                endpoints.push_back(s+(p-s)*d);
//        }
//    }
//    polygon.pop_back();
//
//    for (P p: endpoints) {
//        cout << p << ",";
//    }
//    cout << endpoints[0] << nl;
////    vec<P> tmp{endpoints};
////    endpoints.clear();
////    for (P p: tmp) {
////        for (int i = 0; i < n; i++) {
////
////        }
////    }
//
//    //dijkstra算法求最短长度
//    vector<P> points(polygon.size()+2+endpoints.size());
//    int i=0;
//    points[i++]=g;
//    for(auto p:polygon)points[i++]=p;
//    points[i++]=s;
//    for(auto p:endpoints)points[i++]=p;
//    priority_queue<pair<double,int>, vector<pair<double,int>>, greater<>> que;
//    vector<double> dist(points.size(),1e20);//从g到每个节点的最短路径长度
//    que.emplace(0.0,0);//从0号guard的起始位置出发
//    for(i=que.top().second;i<=n;i=que.top().second){
//        auto d=que.top().first;
//        assert(que.size());
//        que.pop();
//        if(d>=dist[i])continue;
//        dist[i]=d;
//        for(int j=0;j<points.size();++j){
//            if(i==j)continue;
//            auto p1=points[i], p2=points[j];
//            auto ln=(p1-p2).abs();
//            if(ln<eps ||
//               //如果j在poly顶点i关联的两条poly边上,可以直接沿着poly边走到j,不需要后面的包含性判断
//               i>=1&&i<=n && POnLine(polygon[i-1], polygon[(i) % n], p2) ||
//               i>=1&&i<=n && POnLine(polygon[i-1], polygon[(i + n - 2) % n], p2)) {
//                que.emplace(d+ln,j);
//                continue;
//            }
//            //判断p1p2是否被poly完全包含
//            bool fullContained=true;
//            p2=p1+(p2-p1)/ln;
//            for(int k=0;k<polygon.size();++k){
//                auto rd= RayIntersect(p1,p2,polygon[k],polygon[(k+1)%n]);
//                if(rd>eps && rd<ln-eps) {
//                    fullContained=false;
//                    break;
//                }
//            }
//            if(!fullContained)continue;
//            //判断一下p1p2上任意一点是不是在poly内部，来保证p1p2全段都在poly内部。通过p1p2中点画一任意射线, 看穿过poly为奇数还是偶数次
//            int cnt=0;
//            p1 = p1*2/3 + p2/ 3;
//            p2 = p1 + P(cos(10), sin(10));
//            for(int k=0;k<n;++k){
//                auto rd= RayIntersect(p1,p2,polygon[k],polygon[(k+1)%n]);
//                cnt+=(rd>eps);
//            }
//            if(cnt%2==1)
//                que.emplace(d+ln,j);
//        }
//    }
//    setPrec();
//    cout << que.top().first << nl;
//}
#include <bits/stdc++.h>
#define ff first
#define ss second
using namespace std;

typedef long long ll;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;

typedef long double ld;

const int SZ = 110;
const ld eps = 1e-9;

struct point {
    int x, y;
};
point operator+(point a, point b) { return {a.x + b.x, a.y + b.y}; }
point operator-(point a, point b) { return {a.x - b.x, a.y - b.y}; }
int operator*(point a, point b) { return a.x * b.x + a.y * b.y; }
int operator/(point a, point b) { return a.x * b.y - a.y * b.x; }

int ccw(point p, point q, point r) {
    int s = (q - p) / (r - p);
    return (s > 0) - (s < 0);
}
pii cross(point p, point q, point r, point v) {
    if (v / (q - p) < 0 && v / (r - p) < 0) return {-1, 1};
    if (v / (q - p) > 0 && v / (r - p) > 0) return {-1, 2};
    point u = r - q;
    int a = (q - p) / u, b = v / u;
    if (b == 0 || 1ll * a * b < 0) return {-1, 3};
    if (a < 0) a *= -1, b *= -1;
    return {a, b};
}

struct pointd {
    ld x, y;
    pointd() : x(0), y(0) {}
    pointd(point p) : x(p.x), y(p.y) {}
    pointd(ld x, ld y) : x(x), y(y) {}
};
pointd operator+(pointd a, pointd b) { return {a.x + b.x, a.y + b.y}; }
pointd operator-(pointd a, pointd b) { return {a.x - b.x, a.y - b.y}; }
ld operator*(pointd a, pointd b) { return a.x * b.x + a.y * b.y; }
ld operator/(pointd a, pointd b) { return a.x * b.y - a.y * b.x; }
pointd operator*(ld k, pointd a) { return {k * a.x, k * a.y}; }

int ccw(pointd p, pointd q, pointd r) {
    ld s = (q - p) / (r - p);
    return (s > eps) - (s < -eps);
}
bool point_on_line(pointd p, pointd q, pointd r) {
    ld s = (q - p) * (r - p);
    return ccw(p, q, r) == 0 && -eps <= s && s <= (q - p) * (q - p);
}

pointd hypot(pointd p, pointd q, pointd r, bool &flag) {
    ld s = (q - p) * (r - p);
    if (s < -eps || s > (q - p) * (q - p) + eps) {
        flag = false; return r;
    }
    s /= (q - p) * (q - p);

    // cout << "hypot " << "(" << r.x << ", " << r.y << ") to segment (" << p.x << ", " << p.y << ") and (" << q.x << ", " << q.y << ")\n";
    // cout << "(" << (p + s * (q - p)).x << ", " << (p + s * (q - p)).y << ")\n";
    return p + s * (q - p);
}
bool cross(pointd p, pointd q, pointd r, pointd s) {
    if (ccw(p, q, r) == 0 && ccw(p, q, s) == 0) return false;
    if (point_on_line(p, q, r) || point_on_line(p, q, s) || point_on_line(r, s, p) || point_on_line(r, s, q)) return true;
    return ccw(p, q, r) * ccw(p, q, s) < 0 && ccw(r, s, p) * ccw(r, s, q) < 0;
}
ld intersect(pointd p, pointd u, pointd q, pointd v) {
    return ((q - p) / v) / (u / v);
}

int n; point P[SZ];
vector<pair<ld, int>> G[SZ];
ld ds[SZ];

bool polygon_in(pointd p) {
    bool flag = false;
    for (int i = 0; i < n; i++) {
        pointd q = P[i], r = P[(i + 1) % n];
        if (point_on_line(q, r, p)) return true;
    }
    for (int i = 0; i < n; i++) {
        pointd q = P[i], r = P[(i + 1) % n];
        if (q.y < r.y) swap(q, r);
        if (q.y > p.y + eps && r.y > p.y + eps) continue;
        if (q.y < p.y + eps && r.y < p.y + eps) continue;
        if (ccw(q, r, p) <= 0) {
            flag ^= 1;
        }
    }
    return flag;
}

int main() {
    freopen("A.in", "r", stdin);
    cin.tie(0); ios_base::sync_with_stdio(0);
    cout.precision(10);

    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> P[i].x >> P[i].y;
    }
    cin >> P[n].x >> P[n].y;
    cin >> P[n + 1].x >> P[n + 1].y;

    point Q[n];
    for (int i = 0; i < n; i++) Q[i] = P[i];
    auto sgn = [&](point p) {
        return (p.y == 0) ? p.x > 0 : p.y > 0;
    };
    sort(Q, Q + n, [&](point p, point q) {
        if (sgn(p - P[n + 1]) != sgn(q - P[n + 1])) return sgn(p - P[n + 1]);
        else return ccw(P[n + 1], p, q) > 0;
    });

    pointd R[2 * n];
    for (int i = 0; i < n; i++) {
        point v = (Q[i] - P[n + 1]) + (Q[(i + 1) % n] - P[n + 1]);

        // cout << "between " << "(" << Q[i].x << ", " << Q[i].y << ") and (" << Q[(i + 1) % n].x << ", " << Q[(i + 1) % n].y << ")\n";

        pii M = {-1, 0}; int p = -1;
        for (int j = 0; j < n; j++) {
            point q = P[j], r = P[(j + 1) % n];
            pii K = cross(P[n + 1], q, r, v);
            if (K.ff != -1) {
                if (M.ff == -1 || 1ll * M.ff * K.ss > 1ll * M.ss * K.ff) M = K, p = j;
            }
        }
        assert(p != -1);

        pii A = cross(P[n + 1], P[p], P[(p + 1) % n], Q[i] - P[n + 1]);
        pii B = cross(P[n + 1], P[p], P[(p + 1) % n], Q[(i + 1) % n] - P[n + 1]);

        // cout << A.ff << "/" << A.ss << " " << B.ff << "/" << B.ss << '\n';

        R[2 * i] = (ld)A.ff / A.ss * (Q[i] - P[n + 1]) + P[n + 1];
        R[2 * i + 1] = (ld)B.ff / B.ss * (Q[(i + 1) % n] - P[n + 1]) + P[n + 1];


        if (ccw(P[n + 1], R[2 * i], P[n]) > 0 && ccw(P[n + 1], P[n], R[2 * i + 1]) > 0 && ccw(R[2 * i], R[2 * i + 1], P[n]) >= 0) {
            return !(cout << fixed << 0 << '\n');
        }
        if (ccw(P[n + 1], R[2 * i], P[n]) == 0 && (P[n + 1] - P[n]) * (P[n + 1] - P[n]) < (P[n + 1] - R[2 * i]) * (P[n + 1] - R[2 * i]) + eps) {
            return !(cout << fixed << 0 << '\n');
        }
        if (ccw(P[n + 1], R[2 * i + 1], P[n]) == 0 && (P[n + 1] - P[n]) * (P[n + 1] - P[n]) < (P[n + 1] - R[2 * i + 1]) * (P[n + 1] - R[2 * i + 1]) + eps) {
            return !(cout << fixed << 0 << '\n');
        }
    }
    for (int i = 0; i < 2 * n; i++) {
        cout << "(" << R[i].x << ", " << R[i].y << "), ";
    }
    cout << endl;

    vector<pair<pii, pointd>> V;
    for (int i = 0; i <= n; i++) for (int j = 0; j < i; j++) V.push_back({{i, j}, P[j]});
    for (int i = 0; i <= n; i++) for (int j = 0; j < 2 * n; j++) V.push_back({{i, n + 1}, R[j]});
    for (int i = 0; i <= n; i++) for (int j = 0; j < 2 * n; j++) {
            pointd q = R[j], r = R[(j + 1) % (2 * n)];
            if (abs((q - r) * (q - r)) < 1e-12) continue;

            bool flag = true;
            pointd h = hypot(q, r, P[i], flag);
            if (flag) V.push_back({{i, n + 1}, h});
        }

    for (auto [pr, s] : V) {
        auto [i, j] = pr;
        pointd p = P[i];

        vector<ld> A = {0, 1};
        for (int k = 0; k < n; k++) {
            pointd q = P[k], r = P[(k + 1) % n];
            bool f = cross(p, s, q, r);
            if (f) {
                A.push_back(intersect(p, s - p, q, r - q));
            }
        }
        sort(A.begin(), A.end());

        int sz = A.size();
        bool flag = true;
        for (int k = 0; k < sz - 1; k++) {
            ld K = 0.5 * (A[k + 1] + A[k]);
            pointd m = p + K * (s - p);
            if (!polygon_in(m)) {
                flag = false;
                break;
            }
        }
        if (flag) {
            ld d = sqrtl((s - p) * (s - p));
            // cout << "edge between " << i << " and " << j << " -> " << d << '\n';
            G[i].push_back({d, j});
            G[j].push_back({d, i});
        }
    }

    fill(ds, ds + SZ, 1e18);
    priority_queue<pair<ld, int>> pq;
    ds[n] = 0;
    pq.push({0, n});
    while (!pq.empty()) {
        auto [d, v] = pq.top(); pq.pop();
        if (ds[v] != -d) continue;
        for (auto [w, x] : G[v]) {
            if (ds[x] > ds[v] + w) {
                ds[x] = ds[v] + w;
                pq.push({-ds[x], x});
            }
        }
    }
    cout << fixed << ds[n + 1] << '\n';
}
