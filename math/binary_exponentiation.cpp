#include <bits/stdc++.h>
using namespace std;

// O(log(b))
int binpow(int a, int b) {
    int res = 1;
    while (b) {
        if (b & 1)
            res *= a;
        a *= a;
        b >>= 1;
    }
    return res;
}

// O(log(b))
int binpow(int a, int b, const int m) {
    a %= m;
    int res = 1;
    while (b) {
        if (b & 1)
            res = res * a % m;
        a = a * a % m;
        b >>= 1;
    }
    return res;
}