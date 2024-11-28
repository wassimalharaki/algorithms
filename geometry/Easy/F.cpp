#include <cmath>
#include <cassert>
#include <cstdio>

const double PI = acos(-1.0);


double f(double t) {
    return (PI - t) / (1 + cos(t));
}

double getans(double t) {
    return tan(t / 2.0);  // sqrt(2 / (1 + cos(t)) - 1);                                                                                   
}

double argmin() {
    double start = 0, end = PI / 2;
    for (int i = 0; i < 2000; ++i) {
        double mid1 = (end - start) / 3 + start;
        double mid2 = end - (end - start) / 3;
        // if (mid2 - mid1 < 1e-12) break;
        if (f(mid1) < f(mid2)) {
            end = mid2;
        } else {
            start = mid1;
        }
    }
    return (start + end) / 2;
}

int main() {
    int T;
    long long u;
    freopen("help.in", "r", stdin);
    assert(scanf("%d", &T) > 0);
    const double t = argmin();
    const double ans = getans(t);
    const double len = f(t);
    for (int kase = 0; kase < T; ++kase) {
        assert(scanf("%lld", &u) > 0);
        printf("%d %d\n", (int) (u * ans), (int) (u * len));
    }
    return 0;
}