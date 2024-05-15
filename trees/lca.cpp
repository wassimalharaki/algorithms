#include <bits/stdc++.h>
using namespace std;
#define v vector
 
struct LCA {
    using ai2 = array<int, 2>;
    v<int> in;
    v<v<ai2>> d;

    LCA(v<v<int>>& adj, int root = 0) {
        in.resize(adj.size());
 
        v<ai2> path;
        auto dfs = [&](int u, int p, int h, auto&& dfs) -> void {
            in[u] = path.size();
            path.push_back({u, h});
 
            for (int& i : adj[u])
                if (i != p) {
                    dfs(i, u, h + 1, dfs);
                    path.push_back({u, h});
                }
        };
        dfs(root, -1, 0, dfs);
        build(path);
    }

    ai2 op(ai2& l, ai2& r) {
        return l[1] < r[1] ? l : r;
    }

    void build(v<ai2>& a) {
        int n = a.size(), k = 64 - __builtin_clzll(n);
        d = v<v<ai2>>(k, v<ai2>(n));
        copy(a.begin(), a.end(), d[0].begin());
 
        for (int i = 1; i <= k; i++)
            for (int j = 0; j + (1 << i) <= n; j++)
                d[i][j] = op(d[i - 1][j], d[i - 1][j + (1 << (i - 1))]);
    }

    int prod(int a, int b) {
        int l = in[a], r = in[b];
        if (l > r) swap(l, r);
        int i = 63 - __builtin_clzll(++r - l);
        return op(d[i][l], d[i][r - (1 << i)])[0];
    }
};