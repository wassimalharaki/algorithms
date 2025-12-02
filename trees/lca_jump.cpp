// O(n), O(log(n))
struct LCA {
    struct node { int d, jump, p; };
    vector<node> go;

    LCA(const vector<vector<int>>& adj, int root = 0) {
        go.resize(adj.size());
        go[root] = {0, root, root};

        auto dfs = [&](int u, int p, auto&& self) -> void {
            for (const int& i : adj[u])
                if (i != p) {
                    add_leaf(i, u);
                    self(i, u, self);
                }
        };
        dfs(root, -1, dfs);
    }

    LCA(int n, int root = 0) {
        go.resize(n);
        go[root] = {0, root, root};
    }

    void add_leaf(int u, int p) {
        go[u].p = p;
        go[u].d = go[p].d + 1;

        int p2 = go[p].jump;
        if (go[p].d == 2 * go[p2].d - go[go[p2].jump].d)
            go[u].jump = go[p2].jump;
        else
            go[u].jump = p;
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