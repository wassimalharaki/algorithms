#include <bits/stdc++.h>
using namespace std;
#define v vector

// O(V + E)
v<int> get_cycle(v<v<int>>& adj) {
    int n = adj.size(), l = -1, r = -1;
    v<int> p(n, -1);
    v<char> c(n);
    
    auto dfs = [&](int u, auto&& dfs) -> bool {
        c[u] = 1;
        for (int& i : adj[u])
            if (c[i] == 0) {
                p[i] = u;
                if (dfs(i, dfs))
                    return 1;
            }
            else if (c[i] == 1) {
                l = i;
                r = u;
                return 1;
            }
        c[u] = 2;
        return 0;
    };

    for (int i = 0; i < n; i++)
        if (c[i] == 0 and dfs(i, dfs))
            break;
    if (r == -1) return {};

    v<int> cycle{l};
    for (int u = r; u != l; u = p[u])
        cycle.push_back(u);
    cycle.push_back(l);
    reverse(cycle.begin(), cycle.end());

    return cycle;
}