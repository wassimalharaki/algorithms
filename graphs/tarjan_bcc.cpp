// O(V + E)
vector<int> comp, comps;
void tarjan(vector<vector<int>>& adj) {
    int n = adj.size(), curr = 0;
    comps.clear();
    comp.assign(n, -1);
    vector<int> disc(n), vis;

    auto dfs = [&](int u, int p, auto&& dfs) -> int {
        int low = disc[u] = ++curr;
        vis.push_back(u);

        for (int& i : adj[u])
            if (i != p and comp[i] == -1)
                low = min(low, disc[i] ?: dfs(i, u, dfs));

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
        if (!disc[i]) dfs(i, -1, dfs);
    reverse(comps.begin(), comps.end());
}