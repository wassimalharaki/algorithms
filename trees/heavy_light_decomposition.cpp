// O(n), O(log^2(n))
template <class U, class S, S (*op)(S, S), S (*e)()>
struct HLD {
    vector<int> p, d, tree_id, idx, root;
    vector<U> trees;

    HLD(const vector<vector<int>>& adj, const vector<S>& a, int src = 0) {
        int n = adj.size(), curr_root;
        vector<int> cnt(n);
        vector<S> chain;
        p.resize(n);
        d.resize(n);
        tree_id.resize(n);
        idx.resize(n);
        root.resize(n);

        auto dfs1 = [&](int u, int par, auto&& self) -> void {
            p[u] = par, cnt[u]++;
            d[u] = par == -1 ? 0 : d[par] + 1;
            for (const int& i : adj[u])
                if (i != par) {
                    self(i, u, self);
                    cnt[u] += cnt[i];
                }
        };
        dfs1(src, -1, dfs1);

        auto dfs2 = [&](int u, auto&& self) -> void {
            tree_id[u] = trees.size();
            idx[u] = chain.size();
            if (chain.empty())
                curr_root = u;
            chain.push_back(a[u]);
            root[u] = curr_root;

            if (adj[u].size() == 1 and p[u] != -1) {
                trees.push_back(U(chain));
                chain.clear();
            }
            else {
                int imax_cnt = -1;
                for (const int& i : adj[u])
                    if (i != p[u])
                        if (imax_cnt == -1 or cnt[i] > cnt[imax_cnt])
                            imax_cnt = i;
                self(imax_cnt, self);

                for (const int& i : adj[u])
                    if (i != p[u] and i != imax_cnt)
                        self(i, self);
            }
        };
        dfs2(src, dfs2);
    }

    void set(int i, S x) {
        trees[tree_id[i]].set(idx[i], x);
    }

    S get(int i) {
        return trees[tree_id[i]].get(idx[i]);
    }

    S prod(int a, int b) {
        S ans = e();
        while (tree_id[a] != tree_id[b]) {
            if (d[root[a]] < d[root[b]]) swap(a, b);
            ans = op(ans, trees[tree_id[a]].prod(0, idx[a] + 1));
            a = p[root[a]];
        }
        return op(ans, trees[tree_id[a]].prod(
            min(idx[a], idx[b]), max(idx[a], idx[b]) + 1));
    }
};