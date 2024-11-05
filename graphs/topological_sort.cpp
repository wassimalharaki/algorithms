// O(V + E)
vector<int> topsort(const vector<vector<int>>& adj) {
    int n = adj.size();
    vector<char> vis(n);
    vector<int> order;

    auto dfs = [&](int u, auto&& self) -> void {
        vis[u] = 1;
        for (const int& i : adj[u])
            if (not vis[i])
                self(i, self);
        order.push_back(u);
    };

    for (int i = 0; i < n; i++)
        if (not vis[i]) dfs(i, dfs);
    reverse(order.begin(), order.end());

    return order;
}