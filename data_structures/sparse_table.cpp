#include <bits/stdc++.h>
using namespace std;

//O(nlog(n))
template <class S, S (*op)(S, S)>
struct sparse_table {
    const int n = 25;
    vector<vector<S>> d;

    sparse_table(vector<S>& nums) {
        d = vector<vector<S>>(n + 1, vector<S>(nums.size()));
        copy(nums.begin(), nums.end(), d[0].begin());

        for (int i = 1; i <= n; i++)
            for (int j = 0; j + (1 << i) <= nums.size(); j++)
                d[i][j] = op(d[i - 1][j], d[i - 1][j + (1 << (i - 1))]);
    }

    S prod(int l, int r) {
        int i = 63 - __builtin_clzll(r - l);
        return op(d[i][l], d[i][r - (1 << i)]);
    }
};