#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define pii array<int,2>
#define nx(i) (i+1)%n
#define pv(i) (i-1+n)%n

struct Point {
    ll x,y;

    Point operator+(const Point &p) {
        return {x + p.x, y + p.y};
    }

    Point operator-(const Point &p) {
        return {x - p.x, y - p.y};
    }
};

ll cross(Point p1, Point p2) {
    return p1.x * p2.y - p1.y * p2.x;
}

int sign(ll num) {
    if (num < 0) return -1;
    else if (num == 0) return 0;
    else return 1;
}

vector<pii> all_anti_podal(int n, vector<Point> &p) {
    int p1 = 0, p2 = 0; // two "pointers"
    vector<pii> result;

    // parallel edges should't be visited twice
    vector<bool> vis(n, false);

    for (;p1<n;p1++) {
        // the edge that we are going to consider in this iteration
        // the datatype is Point, but it acts as a vector
        Point base = p[nx(p1)] - p[p1];

        // the last condition makes sure that the cross products don't have the same sign
        while (p2 == p1 || p2 == nx(p1) || sign(cross(base, p[nx(p2)] - p[p2])) == sign(cross(base, p[p2] - p[pv(p2)]))) {
            p2 = nx(p2);
        }

        if (vis[p1]) continue;
        vis[p1] = true;

        result.push_back({p1, p2});
        result.push_back({nx(p1), p2});

        // if both edges from p1 and p2 are parallel to each other
        if (cross(base, p[nx(p2)] - p[p2]) == 0) {
            result.push_back({p1, nx(p2)});
            result.push_back({nx(p1), nx(p2)});
            vis[p2] = true;
        }
    }

    return result;
}
