#include <bits/stdc++.h>
using namespace std;
#define v vector

// O(V^2)
bool is_cyclic(v<v<int>>& adj) {
    int n = adj.size();
    
    v<char> color(n);
    auto dfs = [&](int u, auto&& dfs) -> bool {
        color[u] = 1;
        for (int i = 0; i < n; i++) {
            if (not adj[u][i]) continue;
            if ((color[i] == 0 and dfs(i, dfs)) or color[i] == 1)
                return 1;
        }
        color[u] = 2;
        return 0;
    };

    for (int i = 0; i < n; i++)
        if (color[i] == 0 && dfs(i, dfs))
            return 1;
    return 0;
}