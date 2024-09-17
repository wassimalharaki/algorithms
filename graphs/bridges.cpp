vector<array<int, 2>> get_bridges(const vector<vector<int>>& adj) {
    int n = adj.size(), t = 0;
    vector<char> vis(n);
    vector<array<int, 2>> bridges;
    vector<int> tin(n, -1), low(n, -1);

    auto dfs = [&](int u, int p, auto&& dfs) -> void {
        vis[u] = 1;
        tin[u] = low[u] = t++;
        
        bool ps = 0;
        for (const int& i : adj[u])
            if (i == p and not ps)
                ps = 1;
            else if (vis[i])
                low[u] = min(low[u], tin[i]);
            else {
                dfs(i, u, dfs);
                low[u] = min(low[u], low[i]);
                if (low[i] > tin[u])
                    bridges.push_back({u, i});
            }
    };

    for (int i = 0; i < n; ++i)
        if (not vis[i]) dfs(i, -1, dfs);
    return bridges;
}