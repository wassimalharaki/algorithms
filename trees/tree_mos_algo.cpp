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

template <class S, class F>
struct tree_mo {
    struct query { int l, r, i, p; int64_t h; };
    vector<query> q;
    int m = 0;

    int64_t hilbert_order(int l, int r){
        int lg = __lg((r << 1) + 1) | 1;
        int mxN = (1 << lg) - 1;
        int64_t ans = 0;
        for (int i = 1 << (lg - 1); i; i >>= 1) {
            bool rl = l & i, ry = r & i;
            ans = (ans << 2) | (rl ? ry ? 2 : 1 : ry ? 3 : 0);
            if (not rl) {
                if (ry) l ^= mxN, r ^= mxN;
                swap(l, r);
            }
        }
        return ans;
    }

    void prod(int a, int b) {
        q.push_back({a, b, m++, -1, -1});
    }

    vector<F> solve(const vector<vector<int>>& adj, const vector<S>& _a) {
        int n = adj.size();
        LCA lca(adj);
        vector<array<S, 2>> a;
        vector<int> in(n), cnt(n);

        int t = 0;
        auto dfs = [&](int u, int p, auto&& self) -> void {
            in[u] = t++;
            a.push_back({_a[u], u});

            for (const int& i : adj[u])
                if (i != p)
                    self(i, u, self);

            t++;
            a.push_back({_a[u], u});
        };
        dfs(0, -1, dfs);

        for (auto& [l, r, i, p, h] : q) {
            if (in[l] > in[r]) swap(l, r);
            p = in[lca.prod(l, r)];
            l = in[l] + 1, r = in[r] + 1;
            h = hilbert_order(l, r);
        }
        sort(q.begin(), q.end(), [&](auto& x, auto& y) {
            return x.h < y.h;
        });

        auto insert = [&](int i) {

        };

        auto erase = [&](int i) {

        };

        auto handle = [&](int i, bool ins) {
            if (ins) {
                if (++cnt[a[i][1]] == 2)
                    erase(i);
                else
                    insert(i);
            }
            else {
                if (--cnt[a[i][1]] == 1)
                    insert(i);
                else
                    erase(i);
            }
        };

        int l = 0, r = 0;
        vector<F> ans(m);
        for (auto& [_l, _r, i, p, _] : q) {
            while (r < _r) handle(r++, 1);
            while (l > _l) handle(--l, 1);
            while (r > _r) handle(--r, 0);
            while (l < _l) handle(l++, 0);
            handle(p, 1);
            ans[i] = 0;
            handle(p, 0);
        }
        return ans;
    }
};