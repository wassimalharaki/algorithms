#include <bits/stdc++.h>
using namespace std;
#define int long long

signed main() {
    int x = 1234;
    
    //count ones
    __builtin_popcountll(x);

    //index of highest one
    63 - __builtin_clzll(x);

    //index of lowest one
    __builtin_ctzll(x);

    return 0;
}