#include <bits/stdc++.h>
using namespace std;

struct prefix_2d {
    vector<vector<int>> d;

    prefix_2d(vector<vector<int>>& a) {
        int n = a.size();
        d.resize(n + 1, vector<int>(n + 1));

        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
                d[i][j] = d[i - 1][j] + d[i][j - 1]
                    + a[i - 1][j - 1] - d[i - 1][j - 1];
    }

    int prod(int i1, int j1, int i2, int j2) {
        return d[i2][j2] + d[i1][j1]
            - d[i2][j1] - d[i1][j2];
    }
};