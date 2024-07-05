// O(n^2), O(log^2(n))
template <class S, S (*op)(S, S), S (*e)()>
struct segtree_2d {
    int n, m, size_n, size_m, log_n, log_m;
    vector<vector<S>> d;

    segtree_2d(int _n, int _m) :
        segtree_2d(vector(_n, vector<S>(_m, e()))) {}
    
    segtree_2d(const vector<vector<S>>& a) {
        n = a.size();
        m = a[0].size();
        size_n = bit_ceil(n);
        size_m = bit_ceil(m);
        log_n = __builtin_ctz(size_n);
        log_m = __builtin_ctz(size_m);
        d.resize(size_n << 1, vector(size_m << 1, e()));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++)
                d[size_n + i][size_m + j] = a[i][j];
            for (int j = size_m - 1; j >= 1; j--)
                update_inner(size_n + i, j);
        }

        for (int i = size_n - 1; i >= 1; i--)
            for (int j = 0; j < size_m << 1; j++)
                update_outer(i, j);
    }

    void update_inner(int i, int j) {
        d[i][j] = op(d[i][j << 1], d[i][(j << 1) + 1]);
    }

    void update_outer(int i, int j) {
        d[i][j] = op(d[i << 1][j], d[(i << 1) + 1][j]);
    }

    int bit_ceil(int _n) {
        int x = 1;
        while (x < n) x <<= 1;
        return x;
    }

    void set(int r, int c, S x) {
        r += size_n, c += size_m;
        d[r][c] = x;
        for (int i = 1; i <= log_n; i++)
            update_outer(r >> i, c);
        for (int i = 1; i <= log_m; i++) {
            update_inner(r, c >> i);
            for (int j = 1; j <= log_n; j++)
                update_outer(r >> j, c >> i);
        }
    }

    S get(int r, int c) const {
        return d[r + size_n][c + size_m];
    }

    S prod(int i, int l, int r) const {
        S sml = e(), smr = e();
        l += size_m, r += size_m;

        while (l < r) {
            if (l & 1) sml = op(sml, d[i][l++]);
            if (r & 1) smr = op(d[i][--r], smr);
            l >>= 1;
            r >>= 1;
        }
        return op(sml, smr);
    }

    S prod(int r1, int c1, int r2, int c2) const {
        S sml = e(), smr = e();
        r1 += size_n, r2 += size_n;

        while (r1 < r2) {
            if (r1 & 1) sml = op(sml, prod(r1++, c1, c2));
            if (r2 & 1) smr = op(prod(--r2, c1, c2), smr);
            r1 >>= 1;
            r2 >>= 1;
        }
        return op(sml, smr);
    }

    S all_prod() const { return d[1][1]; }
};