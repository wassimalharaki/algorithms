// O(n), O(log(n))
struct DSU {
    int comp;
    vector<int> p;
    vector<array<int, 5>> hist;

    DSU(int n) { p.resize(n, -1); comp = n; }

    int find(int x) {
        return p[x] < 0 ? x : find(p[x]);
    }

    int size(int x) { return - p[find(x)]; }

    bool same_set(int x, int y) {
        return find(x) == find(y);
    }

    bool merge(int x, int y) {
        x = find(x); y = find(y);
        hist.push_back({x, p[x], y, p[y], comp});
        if (x == y) return false;
        if (p[x] > p[y]) swap(x, y);
        p[x] += p[y]; p[y] = x;
        comp--; return true;
    }

    void rollback() {
        auto [x, sx, y, sy, c] = hist.back();
        p[x] = sx; p[y] = sy; comp = c;
        hist.pop_back();
    }
};

// O(qlog^2(q))
struct dc_graph {
    int n, size, t;
    vector<char> q;
    vector<array<int, 3>> e; 
    vector<vector<array<int, 2>>> d;

    dc_graph(int _n, int _q) {
        t = 0, n = _n;
        size = q <= 1 ? 1 : 1 << (1 + __lg(q - 1));;
        q.resize(size);
        d.resize(size << 1);
    }

    void add_edge(int a, int b) {
        if (a > b) swap(a, b);
        e.push_back({a, b, t++});
    }

    void erase_edge(int a, int b) {
        if (a > b) swap(a, b);
        e.push_back({a, b, t++});
    }

    void query() { q[t++] = 1; }

    void apply(int l, int r, array<int, 2> x) {
        l += size, r += size;
        while (l < r) {
            if (l & 1) d[l++].push_back(x);
            if (r & 1) d[--r].push_back(x);
            l >>= 1; r >>= 1;
        }
    }

    vector<int> solve() {
        int m = e.size();
        sort(e.begin(), e.end());

        for (int i = 0; i < m; i++)
            if (i == m - 1)
                apply(e[i][2], t, {e[i][0], e[i][1]});
            else if (e[i][0] == e[i + 1][0] and e[i][1] == e[i + 1][1])
                apply(e[i][2], e[i + 1][2], {e[i][0], e[i][1]}), i++;
            else
                apply(e[i][2], t, {e[i][0], e[i][1]});

        DSU ds(n);
        vector<int> ans;
        auto dfs = [&](int i, auto&& self) -> void {
            for (auto& [a, b] : d[i])
                ds.merge(a, b);

            if (i < size) {
                self(i << 1, self);
                self((i << 1) + 1, self);
            }
            else if (q[i - size])
                ans.push_back(ds.comp);
            
            for (int j = 0; j < (int) d[i].size(); j++)
                ds.rollback();
        };
        dfs(1, dfs);

        return ans;
    }
};