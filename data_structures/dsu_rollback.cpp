// O(n), O(1)
struct DSU {
    vector<int> p;
    vector<array<int, 4>> history;

    DSU(int n) { p.resize(n, -1); }

    int find(int x) {
        return p[x] < 0 ? x : find(p[x]);
    }

    int size(int x) { return - p[find(x)]; }

    bool same_set(int x, int y) {
        return find(x) == find(y);
    }

    bool merge(int x, int y) {
        x = find(x); y = find(y);
        history.push_back({x, p[x], y, p[y]});
        if (x == y) return false;
        if (p[x] > p[y]) swap(x, y);
        p[x] += p[y]; p[y] = x;
        return true;
    }

    void rollback() {
        auto [x, sx, y, sy] = history.back();
        p[x] = sx; p[y] = sy;
        history.pop_back();
    }
};