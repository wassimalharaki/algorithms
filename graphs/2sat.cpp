// O(V + E)
struct two_sat {
    int n;
    vector<bool> ans;
    vector<vector<int>> adj;

    two_sat(int _n) {
        n = _n;
        ans.resize(_n);
        adj.resize(2 * _n);
    }

    void add_clause(int i, bool f, int j, bool g) {
        adj[2 * i + !f].push_back(2 * j + g);
        adj[2 * j + !g].push_back(2 * i + f);
    }

    vector<int> tarjan() {
        int curr = 0, grp_id = 0;
        vector<int> disc(2 * n), id(2 * n, -1), vis;

        auto dfs = [&](int u, auto&& dfs) -> int {
            int low = disc[u] = ++curr;
            vis.push_back(u);

            for (int& i : adj[u])
                if (id[i] == -1)
                    low = min(low, disc[i] ?: dfs(i, dfs));

            if (low == disc[u]) {
                for (int i = -1; i != u;) {
                    i = vis.back();
                    id[i] = grp_id;
                    vis.pop_back();
                }
                grp_id++;
            }
            return low;
        };

        for (int i = 0; i < 2 * n; i++)
            if (!disc[i]) dfs(i, dfs);
        return id;
    }

    bool satisfiable() {
        vector<int> id = tarjan();
        for (int i = 0; i < n; i++) {
            if (id[2 * i] == id[2 * i + 1]) return 0;
            ans[i] = id[2 * i] > id[2 * i + 1];
        }
        return 1;
    }

    vector<bool> answer() { return ans; }
};