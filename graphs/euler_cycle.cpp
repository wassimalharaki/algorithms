#include <bits/stdc++.h>
using namespace std;
#define v vector

v<int> euler_cycle(v<v<array<int, 2>>> adj, int src = 0) {
    int m = 0;
    for (auto& a : adj) {
        if (a.size() & 1)
            return {};
        m += a.size();
    }
    m >>= 1;

    v<char> vis(m);
    v<int> path;
    auto dfs = [&](int u, auto&& dfs) -> void {
        while (adj[u].size()) {
            auto [a, i] = adj[u].back();
            adj[u].pop_back();
            if (not vis[i]) {
                vis[i] = 1;
                dfs(a, dfs);
            }
        }
        path.push_back(u);
    };
    dfs(src, dfs);

    if (path.size() == m + 1)
        return path;
    else
        return {};
}