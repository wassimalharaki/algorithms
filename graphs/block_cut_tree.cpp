struct block_cut_tree {
    vector<int> comp;
    vector<char> is_cutpoint;
    vector<vector<int>> adj;

    block_cut_tree(const vector<vector<int>>& _adj) {
        int n = _adj.size(), curr = 0;
        comp.resize(n);
        is_cutpoint.resize(n);
        vector<vector<int>> comps;
        vector<int> disc(n), low(n), vis;

        auto dfs = [&](int u, int p, auto&& self) -> void {
            disc[u] = low[u] = ++curr;
            vis.push_back(u);

            for (const int& i : _adj[u]) {
                if (i == p) continue;
                if (disc[i]) {
                    low[u] = min(low[u], disc[i]);
                    continue;
                }

                self(i, u, self);
                low[u] = min(low[u], low[i]);
                if (low[i] >= disc[u]) {
                    is_cutpoint[u] = disc[u] > 1 || disc[i] > 2;

                    comps.push_back({u});
                    while (comps.back().back() != i) {
                        comps.back().push_back(vis.back());
                        vis.pop_back();
                    }
                }
            }
        };

        for (int i = 0; i < n; i++, curr = 0)
            if (not disc[i]) dfs(i, -1, dfs);

        int m = 0;
        for (int i = 0; i < n; i++)
            if (is_cutpoint[i])
                comp[i] = m++;

        adj.resize(m + comps.size());
        for (auto& a : comps) {
            for (int& u : a)
                if (not is_cutpoint[u])
                    comp[u] = m;
                else {
                    adj[comp[u]].push_back(m);
                    adj[m].push_back(comp[u]);
                }
            m++;
        }
    }
};