#include <bits/stdc++.h>
using namespace std;
#define v vector

// O(V + E)
v<int> topsort(v<v<int>>& adj) {
    int n = adj.size();
    v<char> vis(n);
    v<int> order;

    auto dfs = [&](int u, auto&& dfs) -> void {
        vis[u] = 1;
        for (int& i : adj[u])
            if (not vis[i])
                dfs(i, dfs);
        order.push_back(u);
    };

    for (int i = 0; i < n; i++)
        if (not vis[i]) dfs(i, dfs);
    reverse(order.begin(), order.end());

    return order;
}