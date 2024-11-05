vector<int> get_cutpoints(const vector<vector<int>>& adj) {
    int n = adj.size(), t = 0;
    vector<char> vis(n), is_cutpoint(n);
    vector<int> tin(n), low(n);

    auto dfs = [&](int u, int p, auto&& self) -> void {
        vis[u] = 1;
        tin[u] = low[u] = t++;

        int cnt = 0;
        for (const int& i : adj[u])
            if (i == p)
                continue;
            else if (vis[i])
                low[u] = min(low[u], tin[i]);
            else {
                self(i, u, self);
                low[u] = min(low[u], low[i]);
                if (low[i] >= tin[u] and p != -1)
                    is_cutpoint[u] = 1;
                cnt++;
            }

        if(p == -1 and cnt > 1)
            is_cutpoint[u] = 1;
    };

    for (int i = 0; i < n; i++)
        if (not vis[i]) dfs(i, -1, dfs);

    vector<int> cutpoints;
    for (int i = 0; i < n; i++)
        if (is_cutpoint[i])
            cutpoints.push_back(i);
    return cutpoints;
}