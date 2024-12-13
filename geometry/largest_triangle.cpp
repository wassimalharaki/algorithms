#include "./point.cpp"

// convex polygons only.
// O(n), does not work for certain class of convex polygons
void o1() {
    // points must be sorted in cw or ccw
    vec<P> p;
    int n = 0;
    auto area = [&](int i, int j, int k) {
        double a = p[i].dist(p[j]);
        double b = p[i].dist(p[j]);
        double c = p[j].dist(p[k]);
        double s = (a+b+c)/2.0;
        double area = sqrt( s*(s-a)*(s-b)*(s-c) );
        return area;
    };

    int A = 0, B = 1, C = 2;
    int bA = A, bB = B, bC = C;

    while (true) { // loop A
        while (true) { // loop B
            while (lse(area(A, B, C), area(A, B, (C + 1) % n))) { // loop C
                C = (C + 1) % n;
            }
            if (lse(area(A, B, C), area(A, (B + 1) % n, C))) B = (B + 1) % n;
            else break;
        }
        if (ls(area(bA, bB, bC), area(A, B, C))) bA = A, bB = B, bC = C;
        A = (A + 1) % n;
        B = (C + 1) % n;
        C = (C + 1) % n;
        if (A == 0) break;
    }
}

// O(n ^ 2) works for any case
void o2() {
    vec<P> p;
    int n = 10;
    atype res = 0;
    auto area = [&](int i, int j, int k) {
        double a = p[i].dist(p[j]);
        double b = p[i].dist(p[j]);
        double c = p[j].dist(p[k]);
        double s = (a+b+c)/2.0;
        double area = sqrt( s*(s-a)*(s-b)*(s-c) );
        return area;
    };

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            int k = (j + 1) % n;
            atype cur = 0;
            while (true) {
                atype a1 = area(i, j, k), a2 = area(i, j, (k + 1) % n);
                if (ls(a1, a2)) {
                    k = (k + 1) % n;
                } else {
                    cur = a1;
                    break;
                }
            }
            res = max(res, cur);
        }
    }
}

// divide and conquer case
