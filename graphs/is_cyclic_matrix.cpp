// O(V^2)
bool is_cyclic(const vector<vector<int>>& adj) {
    int n = adj.size();

    vector<char> c(n);
    auto dfs = [&](int u, auto&& self) -> bool {
        c[u] = 1;
        for (int i = 0; i < n; i++) {
            if (not adj[u][i]) continue;
            if ((c[i] == 0 and self(i, self)) or c[i] == 1)
                return 1;
        }
        c[u] = 2;
        return 0;
    };

    for (int i = 0; i < n; i++)
        if (c[i] == 0 && dfs(i, dfs))
            return 1;
    return 0;
}