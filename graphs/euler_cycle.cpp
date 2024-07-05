// O(E)
vector<int> euler_cycle(vector<vector<array<int, 2>>> adj, int src = 0) {
    int m = 0;
    for (auto& a : adj) {
        if (a.size() & 1)
            return {};
        m += a.size();
    }
    m >>= 1;

    vector<char> vis(m);
    vector<int> path;
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