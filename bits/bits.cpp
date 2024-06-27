#include <bits/stdc++.h>
using namespace std;
#define int long long

signed main() {
    int x = 1234;
    
    //count ones
    __builtin_popcountll(x);

    //most significant bit
    63 - __builtin_clzll(x);
    __lg(x);

    //least significant bit
    __builtin_ctzll(x);

    return 0;
}