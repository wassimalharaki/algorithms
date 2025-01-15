// O(n), O(1)
struct DSU {
    vector<int> p;

    DSU(int n) { p.resize(n, -1); }

    int find(int x) {
        return p[x] < 0 ? x : p[x] = find(p[x]);
    }

    bool merge(int x, int y) {
        x = find(x); y = find(y);
        if (x == y) return false;
        if (p[x] > p[y]) swap(x, y);
        p[x] += p[y]; p[y] = x;
        return true;
    }
};

// O(Elog(E))
int kruskal(int n, vector<array<int, 3>>& edges) {
    sort(edges.begin(), edges.end(), [](auto& x, auto& y) {
        return x[2] < y[2];
    });

    DSU ds(n);
    int cost = 0, j = 0;
    // vector<array<int, 3>> tree(n - 1);
    for (int i = 0; j < n - 1 and i < (int) edges.size(); i++) {
        auto [a, b, c] = edges[i];
        if (not ds.merge(a, b)) continue;
        // tree[j] = {a, b, c};
        cost += edges[i][0], j++;
    }
    return (j == n - 1 ? cost : -1);
}