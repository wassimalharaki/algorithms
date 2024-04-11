#include <bits/stdc++.h>
using namespace std;

//ax + by = gcd(a, b)
int gcd(int a, int b, int& x, int& y) {
    x = 1, y = 0;
    int x1 = 0, y1 = 1;
    while (b) {
        int q = a / b;
        tie(x, x1) = make_pair(x1, x - q * x1);
        tie(y, y1) = make_pair(y1, y - q * y1);
        tie(a, b) = make_pair(b, a - q * b);
    }
    return a; //gcd
}

//ax + by = c
//c must be = gcd(a, b) * k
//CHECK
//a != 0 and b != 0 (either no solution or infinite nb of solutions)
bool find_any_solution(int a, int b, int c, int &x0, int &y0) {
    int g = gcd(abs(a), abs(b), x0, y0);
    if (c % g) return 0;
    c /= g, x0 *= c, y0 *= c;
    if (a < 0) x0 = -x0;
    if (b < 0) y0 = -y0;
    return 1;
}