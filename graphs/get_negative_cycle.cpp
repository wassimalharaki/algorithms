#include <bits/stdc++.h>
using namespace std;
#define v vector

// O(VE)
v<int> get_neg_cycle(int n, v<array<int, 3>>& edges) {
    v<int> d(n), p(n, -1);
    int x;
    for (int i = 0; i < n; i++) {
        x = -1;
        for (array<int, 3>& j : edges)
            if (d[j[0]] + j[2] < d[j[1]]) {
                x = j[1];
                d[x] = d[j[0]] + j[2];
                p[x] = j[0];
            }
    }
    if (x == -1) return {};

    for (int i = 0; i < n; i++, x = p[x]);
    v<int> cycle{x};
    for (int u = p[x]; u != x; u = p[u])
        cycle.push_back(u);
    cycle.push_back(x);
    reverse(cycle.begin(), cycle.end());

    return cycle;
}