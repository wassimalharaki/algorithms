vector<int> centroid_decomp(const vector<vector<int>>& adj) {
    vector<int> sz(adj.size()), par(adj.size()), dead(adj.size());

    auto find_size = [&](int u, int p, auto&& self) -> int {
        sz[u] = 1;
        for (const int& i : adj[u])
            if (i != p and not dead[i])
                sz[u] += self(i, u, self);
        return sz[u];
    };

    auto find_centroid = [&](int u, int p, int n, auto&& self) -> int {
        for (const int& i : adj[u])
            if (i != p and not dead[i] and sz[i] > n / 2)
                return self(i, u, n, self);
        return u;
    };

    auto build = [&](int u, int p, auto&& self) -> void {
        int n = find_size(u, p, find_size);
        int cent = find_centroid(u, p, n, find_centroid);
        par[cent] = p;
        dead[cent] = 1;

        for (const int& i : adj[cent])
            if (i != p and not dead[i])
                self(i, cent, self);
    };
    build(0, -1, build);
    return par;
}