#include <bits/stdc++.h>
using namespace std;

#define pb push_back
#define mp make_pair
#define fs first
#define sc second
#define double long double

const double pi = acos(-1.0);
const double eps1 = 1e-7;
const double eps2 = 1e-14;

double sqr(double x) {
    return x * x;
}

struct pt {
    double x, y;

    pt(double x_ = 0.0, double y_ = 0.0) : x(x_), y(y_) {}
};

pt operator + (const pt& a, const pt& b) {
    return pt(a.x + b.x, a.y + b.y);
}

pt operator - (const pt& a, const pt& b) {
    return pt(a.x - b.x, a.y - b.y);
}

pt operator * (const pt& a, double k) {
    return pt(a.x * k, a.y * k);
}

struct line {
    double a, b, c;
    line(double a_ = 1.0, double b_ = 0.0, double c_ = 0.0): a(a_), b(b_), c(c_) {
        norm();
    }

    line(const pt& p1, const pt& p2) {
        a = p2.y - p1.y;
        b = p1.x - p2.x;
        c = -a * p1.x - b * p1.y;

        norm();
    }

    void norm() {
        double k = sqrt(sqr(a) + sqr(b));
        a /= k;
        b /= k;
        c /= k;
    }

    double dist(const pt& p) const {
        return a * p.x + b * p.y + c;
    }
};
bool cmp (pt a, pt b) {
    return (a.x < b.x) || (a.x == b.x && a.y < b.y);
}

bool cw (pt a, pt b, pt c) {
    return a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y) < -eps2;
}

bool ccw (pt a, pt b, pt c) {
    return a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y) > eps2;
}

void convex_hull (vector<pt> & a) {
    if (a.size() == 1)  return;
    if (a.size() == 2) return;
    sort (a.begin(), a.end(), &cmp);
    pt p1 = a[0],  p2 = a.back();
    vector<pt> up, down;
    up.push_back (p1);
    down.push_back (p1);
    for (size_t i=1; i<a.size(); ++i) {
        if (i==a.size()-1 || cw (p1, a[i], p2)) {
            while (up.size()>=2 && !cw (up[up.size()-2], up[up.size()-1], a[i]))
                up.pop_back();
            up.push_back (a[i]);
        }
        if (i==a.size()-1 || ccw (p1, a[i], p2)) {
            while (down.size()>=2 && !ccw (down[down.size()-2], down[down.size()-1], a[i]))
                down.pop_back();
            down.push_back (a[i]);
        }
    }
    a.clear();
    for (size_t i=0; i<up.size(); ++i)
        a.push_back (up[i]);
    for (size_t i=down.size()-2; i>0; --i)
        a.push_back (down[i]);
}

int sign(double x) {
    if (fabs(x) < eps2)
        return 0;
    if (x > 0)
        return 1;
    return -1;
}

int main() {
    //freopen("input.txt", "r", stdin);
    //freopen("output.txt", "w", stdout);

    int n, m;

    cin >> n >> m;

    vector <pt> red(n), blue(m);
    for (int i = 0; i < n; i++) {
        cin >> red[i].x >> red[i].y;
    }

    for (int i = 0; i < m; i++) {
        cin >> blue[i].x >> blue[i].y;
    }

    if (m == 1) {
        cout << -1 << endl;
        return 0;
    }
    vector <pt> cvx_blue = blue;
    convex_hull(cvx_blue);
    int k = cvx_blue.size();

    for (int i = 0; i < n; i++) {
        bool flag = true;
        for (int j = 0; j < k; j++) {
            line hl(cvx_blue[j % k], cvx_blue[(j + 1) % k]);
            if (hl.dist(red[i]) > eps1) {
                flag = false;
            }
        }

        if (flag) {
            cout << -1 << endl;

            return 0;
        }

    }


    double ans = 0.0;

    for (int i = 0; i < m; i++) {
        vector <pt> new_blue;
        for (int j = 0; j < i; j++) {
            new_blue.pb(blue[j] - blue[i]);
            double d = sqr(new_blue.back().x) + sqr(new_blue.back().y);
            new_blue.back().x /= d;
            new_blue.back().y /= d;
        }

        new_blue.pb(pt(0.0, 0.0));

        for (int j = i + 1; j < m; j++) {
            new_blue.pb(blue[j] - blue[i]);
            double d = sqr(new_blue.back().x) + sqr(new_blue.back().y);
            new_blue.back().x /= d;
            new_blue.back().y /= d;
        }

        vector <pt> new_red;
        for (int j = 0; j < n; j++) {
            new_red.pb(red[j] - blue[i]);
            double d = sqr(new_red.back().x) + sqr(new_red.back().y);
            new_red.back().x /= d;
            new_red.back().y /= d;
        }

        convex_hull(new_blue);

        int k = new_blue.size();
        int p = new_red.size();

        for (int j = 0; j < k; j++) {
            bool flag = false;
            line hl(new_blue[j], new_blue[(j + 1) % k]);

            for (int q = 0; q < p; q++) {
                if (hl.dist(new_red[q]) < eps2)
                    flag = true;
            }

            if (flag) {
                if (abs(hl.dist(pt(0.0, 0.0))) < eps2) {
                    cout << -1 << endl;

                    return 0;
                }
                //            cerr << hl.dist(pt(0.0, 0.0)) << endl;
                /*
                if (fabs(hl.dist(pt(0.0, 0.0))) < eps) {
                    cout << -1 << endl;
                    return 0;
                }
                */
                ans = max(ans, 1.0 / abs(hl.dist(pt(0.0, 0.0))));
            }
        }

        for (int q = 0; q < p; q++) {
            for (int j = 0; j < k; j++) {
                line hl(new_red[q], new_blue[j]);
                if (sign(hl.dist(new_blue[(j + 1) % k])) * sign(hl.dist(new_blue[(j + k - 1) % k])) >= 0) {
                    if (abs(hl.dist(pt(0.0, 0.0))) < eps2) {
                        cout << -1 << endl;

                        return 0;
                    }

                    //                  cerr << hl.dist(pt(0.0, 0.0)) << endl;
                    ans = max(ans, 1.0 / abs(hl.dist(pt(0.0, 0.0))));
                }
            }
        }
    }

    cout.precision(20);
    cout << ans / 2.0 << endl;

    return 0;
}