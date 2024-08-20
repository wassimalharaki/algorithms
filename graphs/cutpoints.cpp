void get_cutpoints(vector<vector<int>>& adj) {
    int n = adj.size(), t = 0;
    vector<char> vis;
    vector<int> tin, low, cutpoints;

    auto dfs = [&](int u, int p, auto&& dfs) -> void {
        vis[u] = 1;
        tin[u] = low[u] = t++;

        int cnt = 0;
        for (int i : adj[u])
            if (i == p)
                continue;
            else if (vis[i])
                low[u] = min(low[u], tin[i]);
            else {
                dfs(i, u, dfs);
                low[u] = min(low[u], low[i]);
                if (low[i] >= tin[u] and p != -1)
                    cutpoints.push_back(u);
                cnt++;
            }

        if(p == -1 && cnt > 1)
            cutpoints.push_back(u);
    };

    for (int i = 0; i < n; ++i)
        if (!vis[i]) dfs(i, -1, dfs);
}