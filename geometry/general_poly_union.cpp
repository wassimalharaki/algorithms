#include<bits/stdc++.h>
using namespace std;

typedef long double ld;
const ld eps = 1e-8;
int dcmp(ld x) {
    if(abs(x) < eps) return 0;
    else return x < 0 ? -1 : 1;
}
struct Pt {
    ld x, y;
    Pt(ld _x=0, ld _y=0):x(_x), y(_y) {}
    Pt operator+(const Pt &a) const {
        return Pt(x+a.x, y+a.y);  }
    Pt operator-(const Pt &a) const {
        return Pt(x-a.x, y-a.y);  }
    Pt operator*(const ld &a) const {
        return Pt(x*a, y*a);  }
    Pt operator/(const ld &a) const {
        return Pt(x/a, y/a);  }
    ld operator*(const Pt &a) const {
        return x*a.x + y*a.y;  }
    ld operator^(const Pt &a) const {
        return x*a.y - y*a.x;  }
    bool operator<(const Pt &a) const {
        return dcmp(x-a.x) < 0 || (dcmp(x-a.x) == 0 && dcmp(y-a.y) < 0); }
    bool operator>(const Pt &a) const {
        return dcmp(x-a.x) > 0 || (dcmp(x-a.x) == 0 && dcmp(y-a.y) > 0); }
    bool operator==(const Pt &a) const {
        return dcmp(x-a.x) == 0 && dcmp(y-a.y) == 0;  }
};
const Pt O = {0, 0};
std::ostream& operator<<(std::ostream& os, const Pt& pt) {return os << "(" << pt.x << ", " << pt.y << ")";}

struct Line {
    Pt s, e, v; // start, end, end-start
    ld ang;
    Line(Pt _s=Pt(0, 0), Pt _e=Pt(0, 0)):s(_s), e(_e) { v = e-s; ang = atan2(v.y, v.x); }
    bool operator<(const Line &L) const {
        return ang < L.ang;
    } };

#define all(a) a.begin(),a.end()
int sgn(ld x) { return (x > -eps) - (x < eps); }
ld abs(Pt a) { return sqrt(a * a); }
ld abs2(Pt a) { return a * a; }
ld ori(Pt a, Pt b, Pt c) { return (b - a) ^ (c - a); }
ld arg(Pt x) { return atan2(x.y, x.x); }
bool argcmp(const Pt &a, const Pt &b) { // arg(a) < arg(b)
    int f = (Pt{a.y, -a.x} > Pt{} ? 1 : -1) * (a != Pt{});
    int g = (Pt{b.y, -b.x} > Pt{} ? 1 : -1) * (b != Pt{});
    return f == g ? (a ^ b) > 0 : f < g;
}
int PtSide(Pt p, Line L) {
    return sgn(ori(L.s, L.e, p));
}
Pt LineInter(Line l, Line m) {
    ld s = ori(m.s, m.e, l.s), t = ori(m.s, m.e, l.e);
    return (l.e * s - l.s * t) / (s - t);
}

// Area[i] : area covered by at least i polygon
// P should be counterclockwise
vector<ld> PolyUnion(const vector<vector<Pt>> &P) {
    const int n = P.size();
    vector<ld> Area(n + 1);
    vector<Line> Ls;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < P[i].size(); j++)
            Ls.push_back({P[i][j], P[i][(j + 1) % P[i].size()]});
    auto cmp = [&](Line &l, Line &r) {
        Pt u = l.e - l.s, v = r.e - r.s;
        if (argcmp(u, v)) return true;
        if (argcmp(v, u)) return false;
        return PtSide(l.s, r) < 0;
    };
    sort(all(Ls), cmp);
    for (int l = 0, r = 0; l < Ls.size(); l = r) {
        while (r < Ls.size() and !cmp(Ls[l], Ls[r])) r++;
        Line L = Ls[l];
        vector<pair<Pt, int>> event;
        for (auto it : Ls) {
            Pt c = it.s, d = it.e;
            if (sgn((L.s - L.e) ^ (c - d)) != 0) {
                int s1 = PtSide(c, L) == 1;
                int s2 = PtSide(d, L) == 1;
                if (s1 ^ s2) event.emplace_back(LineInter(L, {c, d}), s1 ? 1 : -1);
            } else if (PtSide(c, L) == 0 and sgn((L.s - L.e) * (c - d)) > 0) {
                event.emplace_back(c, 2);
                event.emplace_back(d, -2);
            }
        }
        sort(all(event), [&](auto i, auto j) {
            return (L.s - i.first) * (L.s - L.e) < (L.s - j.first) * (L.s - L.e);
        });
        int cov = 0, tag = 0;
        Pt lst{0, 0};
        for (auto [p, s] : event) {
            if (cov >= tag) {
                Area[cov] += lst ^ p;
                Area[cov - tag] -= lst ^ p;
            }
            if (abs(s) == 1) cov += s;
            else tag += s / 2;
            lst = p;
        }
    }
    for (int i = n - 1; i >= 0; i--) Area[i] += Area[i + 1];
    for (int i = 1; i <= n; i++) Area[i] /= 2;
    return Area;
};
