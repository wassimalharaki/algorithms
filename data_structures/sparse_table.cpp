#include <bits/stdc++.h>
using namespace std;
#define v vector

// O(nlog(n)), O(1)
template <class S, S (*op)(S, S)>
struct sparse_table {
    v<v<S>> d;

    sparse_table(v<S>& nums) {
        int n = nums.size(), k = 64 - __builtin_clzll(n);
        d = v<v<S>>(k, v<S>(n));
        copy(nums.begin(), nums.end(), d[0].begin());

        for (int i = 1; i <= k; i++)
            for (int j = 0; j + (1 << i) <= n; j++)
                d[i][j] = op(d[i - 1][j], d[i - 1][j + (1 << (i - 1))]);
    }

    S prod(int l, int r) {
        int i = 63 - __builtin_clzll(r - l);
        return op(d[i][l], d[i][r - (1 << i)]);
    }
};