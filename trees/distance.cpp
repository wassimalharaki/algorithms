// O(nlog(n)), O(1)
struct LCA {
    using ai2 = array<int, 2>;
    vector<int> in, dep;
    vector<vector<ai2>> d;

    LCA(const vector<vector<int>>& adj, int root = 0) {
        in.resize(adj.size());
        dep.resize(adj.size());

        vector<ai2> path;
        auto dfs = [&](int u, int p, auto&& self) -> void {
            in[u] = path.size();
            dep[u] = p == -1 ? 0 : dep[p] + 1;
            path.push_back({u, dep[u]});

            for (const int& i : adj[u])
                if (i != p) {
                    self(i, u, self);
                    path.push_back({u, dep[u]});
                }
        };
        dfs(root, -1, dfs);
        build(path);
    }

    ai2 op(ai2& l, ai2& r) {
        return l[1] < r[1] ? l : r;
    }

    void build(vector<ai2>& a) {
        int n = a.size(), k = 1 + (n ? __lg(n) : 0);
        d.resize(k, vector<ai2>(n));
        copy(a.begin(), a.end(), d[0].begin());

        for (int i = 1; i <= k; i++)
            for (int j = 0; j + (1 << i) <= n; j++)
                d[i][j] = op(d[i - 1][j], d[i - 1][j + (1 << (i - 1))]);
    }

    int prod(int a, int b) {
        int l = in[a], r = in[b];
        if (l > r) swap(l, r);
        int i = __lg(++r - l);
        return op(d[i][l], d[i][r - (1 << i)])[0];
    }

    int dist(int a, int b) {
        return dep[a] + dep[b] - 2 * dep[prod(a, b)];
    }
};