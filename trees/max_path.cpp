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

// O(nlog(n)), O(1)
struct max_path {
    vector<int> in;
    vector<vector<int>> d;

    max_path(vector<array<int, 3>>& e) {
        int n = e.size() + 1;
        sort(e.begin(), e.end());

        vector<list<int>> lt(n);
        for (int i = 0; i < n; i++)
            lt[i].push_back(i);

        DSU ds(n);
        for (auto& [c, a, b] : e) {
            int pa = ds.find(a), pb = ds.find(b);
            if (ds.p[pa] > ds.p[pb]) swap(pa, pb);
            lt[pa].push_back(c);
            lt[pa].splice(lt[pa].end(), lt[pb]);
            ds.merge(pa, pb);
        }

        int k = 0;
        bool par = 0;
        in.resize(n);
        vector<int> a;
        for (int& x : lt[ds.find(0)]) {
            if (par) a.push_back(x);
            else in[x] = k++;
            par = !par;
        }
        build(a);
    }

    void build(const vector<int>& a) {
        int n = a.size(), k = 1 + (n ? __lg(n) : 0);
        d.resize(k, vector<int>(n));
        copy(a.begin(), a.end(), d[0].begin());

        for (int i = 1; i <= k; i++)
            for (int j = 0; j + (1 << i) <= n; j++)
                d[i][j] = max(d[i - 1][j], d[i - 1][j + (1 << (i - 1))]);
    }

    int prod(int a, int b) {
        if (a == b) return 0;
        int l = in[a], r = in[b];
        if (l > r) swap(l, r);
        int i = __lg(r - l);
        return max(d[i][l], d[i][r - (1 << i)]);
    }
};