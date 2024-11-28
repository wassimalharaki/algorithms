#include <bits/stdc++.h>
using namespace std;
#define fi first
#define se second
template<typename T> using MinHeap = priority_queue<T, vector<T>, greater<T>>;
template<typename T> using MaxHeap = priority_queue<T, vector<T>, less<T>>;
constexpr double eps = 1e-9;
int dSgn(double x) { return x < -eps ? -1 : x > eps ? 1 : 0; }
int dCmp(double x, double y) { return dSgn(x - y);  }
template<int CloseR = 1, int CloseL = 1> int dIn(double x, double l, double r)
{ return (CloseL ? l - eps : l + eps) < x && x < (CloseR ? r + eps : r - eps); }
typedef struct Pnt { double x, y; } Vec;
using cPnt = const Pnt; using cVec = const Vec;
inline Vec operator + (cVec a, cVec b) { return {a.x + b.x, a.y + b.y}; }
inline Vec operator - (cVec a, cVec b) { return {a.x - b.x, a.y - b.y}; }
inline Vec operator * (cVec a, const double b) { return {a.x * b, a.y * b}; }
inline Vec operator / (cVec a, const double b) { return {a.x / b, a.y / b}; }
inline int operator == (cPnt a, cPnt b) { return !dCmp(a.x, b.x) && !dCmp(a.y, b.y); }
inline int operator != (cPnt a, cPnt b) { return dCmp(a.x, b.x) || dCmp(a.y, b.y); }
inline double Dot(cVec a, cVec b) { return a.x * b.x + a.y * b.y; }
inline double Cro(cVec a, cVec b) { return a.x * b.y - a.y * b.x; }
inline double Len2(cVec a) { return Dot(a, a); }
inline double Len(cVec a) { return sqrt(Len2(a)); }
inline istream &operator >> (istream &in, Pnt &p) { return in >> p.x >> p.y; }
inline ostream &operator << (ostream &ou, Pnt &p) { return ou << p.x << " " << p.y; }
template<int Strict = 0> inline int In(cPnt p, cPnt a, cPnt b) {
    if(Strict && (p == a || p == b)) return 0;
    cVec PA = a - p, PB = b - p;
    return !dSgn(Cro(PA, PB)) && dSgn(Dot(PA, PB)) <= 0;
}
struct Polygon : vector<Pnt> {
    Polygon(int n) : vector<Pnt>(n) {}
    int nxtp(int x) const { return x + 1 == (int)size() ? 0 : x + 1; }
    void getXs(vector<Pnt> &ret, cPnt u, cPnt v) const {
        for(int i = 0; i < (int)size(); ++i) {
            cPnt p = (*this)[i], q = (*this)[nxtp(i)];
            cVec UV = v - u, PU = u - p, PV = v - p, PQ = q - p;
            if(!dSgn(Cro(UV, PQ))) continue;
            const double area0 = Cro(PU, PQ), area1 = Cro(PV, PQ);
            cPnt C = u + UV * (area0 / (area0 - area1));
            if(In<1>(C, p, q)) ret.push_back(C);
        }
    }
    int containLine(cPnt u, cPnt v) const {
        static vector<Pnt> buf; buf.clear(); getXs(buf, u, v);
        for(cPnt p : buf) if(In<1>(p, u, v)) return 0;
        double fl = 0, fr = 1; const double len = Len(v - u);
        const auto maintain = [&](cPnt p, cPnt q) {
            const auto &calc = [&](cPnt p) { return Len(p - u) / len; };
            double x = calc(p), y = calc(q); if(x > y) ::swap(x, y);
            if(dCmp(x, fl) != 1) fl = max(fl, y); else fr = min(fr, x);
        };
        for(int i = 0; i < (int)size(); ++i) {
            cPnt p = (*this)[i], z = (*this)[nxtp(i)];
            if(In<1>(p, u, v) || !In<1>(z, u, v)) continue;
            for(int j = nxtp(i);; j = nxtp(j)) {
                cPnt q = (*this)[nxtp(j)]; if(In<1>(q, u, v)) continue;
                const int x = dSgn(Cro(p - u, v - u));
                const int y = dSgn(Cro(q - u, v - u)); if(x == -y) return !x;
                maintain(x ? z : p, y ? (*this)[j] : q); break;
            }
        }
        const double f = (fl + fr) * .5;
        cPnt m = u * (1 - f) + v * f; int result = 0;
        for(int i = 0; i < (int)size(); ++i) {
            cPnt p = (*this)[i], q = (*this)[nxtp(i)];
            if(In(m, p, q)) return 1;
            if(p.y < q.y && dIn<0>(m.y, p.y, q.y))
                result += dSgn(Cro(q - p, m - p)) == 1;
            if(p.y > q.y && dIn<0>(m.y, q.y, p.y))
                result += dSgn(Cro(p - q, m - q)) == 1;
        }
        return result & 1;
    }
};
int main() {
    cin.tie(0)->sync_with_stdio(0);
    int n; cin >> n; Polygon pg(n); Pnt goal;
    for(int i = 0; i < n; ++i) cin >> pg[i];
    vector<pair<Pnt, double>> keys(n + 1);
    cin >> keys[0].fi >> goal, keys[0].se = 1e18;
    for(int i = 1; i <= n; ++i) keys[i] = {pg[i - 1], 1e18};
    const auto &check = [&](cPnt x) {
        int l = 0, r = 0;
        for(int i = 0; i < (int)pg.size(); ++i) {
            cPnt p = pg[i], q = pg[pg.nxtp(i)];
            if(In(p, x, goal) && x != p) {
                const int w = dSgn(Cro(x - goal, q - goal));
                if(w) { w == 1 ? ++l : ++r; if(l && r) return 0; }
            }
            if(In(q, x, goal) && x != q) {
                const int w = dSgn(Cro(x - goal, p - goal));
                if(w) { w == 1 ? ++l : ++r; if(l && r) return 0; }
            }
        }
        return 1;
    };
    const auto &checkpnt = [&](cPnt p) {
        return pg.containLine(p, goal) && check(p); };
    for(int i = 0; i <= n; ++i) if(checkpnt(keys[i].fi)) {
            keys[i].se = 0; if(keys[i].fi == goal) continue;
            cPnt P = keys[i].fi; cVec OP = P - goal;
            static vector<Pnt> buf; pair<double, int> mn = {1e18, -1};
            buf.clear(), pg.getXs(buf, P, goal);
            for(int j = 0; j < (int)buf.size(); ++j)
                if(dSgn(Dot(buf[j] - goal, OP)) == 1)
                    mn = min(mn, {Len2(keys[i].fi - buf[j]), j});
            if(mn.se != -1 && check(buf[mn.se])) keys.emplace_back(buf[mn.se], 0);
        }
    const int z = (int)keys.size();
    for(int i = 0; i < z; ++i) if(keys[i].se == 0) {
            cPnt P = keys[i].fi; cVec OP = P - goal;
            const double len = Len(OP), len2 = Len2(OP);
            for(int j = 0; j < z; ++j) if(keys[j].se != 0) {
                    cPnt Q = keys[j].fi; cVec OQ = Q - goal;
                    const double fr = Dot(OP, OQ) / len2;
                    if(!dIn(fr, 0, 1)) continue;
                    if(pg.containLine(Q, goal + OP * fr))
                        keys[j].se = min(keys[j].se, abs(Cro(OP, OQ)) / len);
                }
        }
    MinHeap<pair<double, int>> Q;
    for(int i = 0; i < z; ++i) Q.emplace(keys[i].se, i);
    while(!Q.empty() && Q.top().se) {
        const auto [d, u] = Q.top(); Q.pop();
        if(keys[u].se != d) continue;
        for(int v = 0; v < z; ++v) {
            const double len = Len(keys[u].fi - keys[v].fi);
            if(d + len < keys[v].se && pg.containLine( \
                keys[u].fi, keys[v].fi)) Q.emplace(keys[v].se = d + len, v);
        }
    }
    printf("%.13lf", keys[0].se);
}