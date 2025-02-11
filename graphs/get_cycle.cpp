// O(V + E)
vector<int> get_cycle(const vector<vector<int>>& adj) {
    int n = adj.size(), l = -1, r = -1;
    vector<int> p(n, -1);
    vector<char> c(n);
    
    auto dfs = [&](int u, auto&& self) -> bool {
        c[u] = 1;
        for (const int& i : adj[u])
            if (c[i] == 0) {
                p[i] = u;
                if (self(i, self))
                    return 1;
            }
            else if (c[i] == 1) {
                l = i;
                r = u;
                return 1;
            }
        c[u] = 2;
        return 0;
    };

    for (int i = 0; i < n; i++)
        if (c[i] == 0 and dfs(i, dfs))
            break;
    if (r == -1) return {};

    vector<int> cycle{l};
    for (int u = r; u != l; u = p[u])
        cycle.push_back(u);
    cycle.push_back(l);
    reverse(cycle.begin(), cycle.end());

    return cycle;
}