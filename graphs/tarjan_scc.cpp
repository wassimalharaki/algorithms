// O(V + E)
vector<int> comp, comps;
void tarjan(const vector<vector<int>>& adj) {
    int n = adj.size(), curr = 0;
    comps.clear();
    comp.assign(n, -1);
    vector<int> disc(n), vis;

    auto dfs = [&](int u, auto&& dfs) -> int {
        int low = disc[u] = ++curr;
        vis.push_back(u);

        for (const int& i : adj[u])
            if (comp[i] == -1)
                low = min(low, disc[i] ?: dfs(i, dfs));

        if (low == disc[u]) {
            comps.push_back(u);
            for (int i = -1; i != u;) {
                i = vis.back();
                comp[i] = u;
                vis.pop_back();
            }
        }
        return low;
    };

    for (int i = 0; i < n; i++)
        if (not disc[i]) dfs(i, dfs);
    reverse(comps.begin(), comps.end());
}