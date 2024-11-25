struct LCA {
    struct node { int d, jump, p; };
    vector<node> go;

    LCA(const vector<vector<int>>& adj, int root = 0) {
        go.resize(adj.size());
        go[root] = {0, root, root};

        auto dfs = [&](int u, int p, auto&& self) -> void {
            for (const int& i : adj[u]) {
                if (i == p) continue;
                go[i].p = u;
                go[i].d = go[u].d + 1;

                int p2 = go[u].jump;
                if (go[u].d == 2 * go[p2].d - go[go[p2].jump].d)
                    go[i].jump = go[p2].jump;
                else
                    go[i].jump = u;
                self(i, u, self);
            }
        };
        dfs(root, -1, dfs);
    }

    int jump(int a, int k) {
        while (k and go[a].d)
            if (go[a].d - k <= go[go[a].jump].d) {
                k -= go[a].d - go[go[a].jump].d;
                a = go[a].jump;
            }
            else {
                k--;
                a = go[a].p;
            }
        return a;
    }

    int prod(int a, int b) {
        if (go[a].d > go[b].d)
            swap(a, b);
        b = jump(b, go[b].d - go[a].d);

        while (a != b)
            if (go[a].jump == go[b].jump)
                a = go[a].p, b = go[b].p;
            else
                a = go[a].jump, b = go[b].jump;
        return a;
    }
};