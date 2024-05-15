#include <bits/stdc++.h>
using namespace std;

// O(VE)
vector<int> bellman_ford(int src, int n, vector<vector<int>>& edges) {
    vector<int> d(n, INT_MAX);
    d[src] = 0;

    for (int i = 0; i < n - 1; i++)
        for (auto& j : edges)
            if (d[j[0]] + j[2] < d[j[1]])
                d[j[1]] = d[j[0]] + j[2];

    // TO DETECT NEGATIVE WEIGHT CYCLES
    for (auto j : edges)
        if (d[j[0]] + j[2] < d[j[1]])
            return {};

    return d;
}