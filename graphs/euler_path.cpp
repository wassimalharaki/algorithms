#include <bits/stdc++.h>
using namespace std;
#define v vector

v<int> euler_path(v<v<int>> adj, int src, int dst) {
    int n = adj.size(), m = 0;

    v<int> in(n), out(n);
    for (int i = 0; i < n; i++)
        for (auto& u : adj[i])
            in[u]++, out[i]++, m++;
    
    if (out[src] != in[src] + 1 or out[dst] != in[dst] - 1)
        return {};
    
    for (int i = 0; i < n; i++)
        if (i != src and i != dst and in[i] != out[i])
            return {};

    v<int> path;
    auto dfs = [&](int u, auto&& dfs) -> void {
        while (adj[u].size()) {
            int i = adj[u].back();
            adj[u].pop_back();
            dfs(i, dfs);
        }
        path.push_back(u);
    };
    dfs(src, dfs);

    reverse(path.begin(), path.end());
    if (path.size() == m + 1 and path.back() == dst)
        return path;
    else
        return {};
}