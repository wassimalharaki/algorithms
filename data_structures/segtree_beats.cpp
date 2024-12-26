// O(n), O(log^2(n))
template <class S>
struct segtree_beats {
    static constexpr S INF = numeric_limits<S>::max() / 2.1;
    struct alignas(32) node {
        S sum = 0, g1 = 0, l1 = 0;
        S g2 = -INF, gc = 1, l2 = INF, lc = 1, add = 0;
    };
    vector<node> d;
    int size, log;

    segtree_beats(const vector<S>& a) {
        int n = a.size();
        size = n <= 1 ? 1 : 1 << (1 + __lg(n - 1));
        log = __builtin_ctz(size);
        d.resize(size << 1);
        for (int i = 0; i < n; i++)
            d[size + i].sum = d[size + i].g1 = d[size + i].l1 = a[i];
        for (int i = size - 1; i >= 1; i--) update(i);
    }

    void range_chmin(int l, int r, S x) { inner_apply<1>(l, r, x); }
    void range_chmax(int l, int r, S x) { inner_apply<2>(l, r, x); }
    void range_add(int l, int r, S x) { inner_apply<3>(l, r, x); }
    void range_update(int l, int r, S x) { inner_apply<4>(l, r, x); }
    S range_min(int l, int r) { return inner_fold<1>(l, r); }
    S range_max(int l, int r) { return inner_fold<2>(l, r); }
    S range_sum(int l, int r) { return inner_fold<3>(l, r); }

    void update(int k) {
        node& p = d[k];
        node& l = d[k * 2 + 0];
        node& r = d[k * 2 + 1];

        p.sum = l.sum + r.sum;

        if (l.g1 == r.g1) {
            p.g1 = l.g1;
            p.g2 = max(l.g2, r.g2);
            p.gc = l.gc + r.gc;
        } else {
            bool f = l.g1 > r.g1;
            p.g1 = f ? l.g1 : r.g1;
            p.gc = f ? l.gc : r.gc;
            p.g2 = max(f ? r.g1 : l.g1, f ? l.g2 : r.g2);
        }

        if (l.l1 == r.l1) {
            p.l1 = l.l1;
            p.l2 = min(l.l2, r.l2);
            p.lc = l.lc + r.lc;
        } else {
            bool f = l.l1 < r.l1;
            p.l1 = f ? l.l1 : r.l1;
            p.lc = f ? l.lc : r.lc;
            p.l2 = min(f ? r.l1 : l.l1, f ? l.l2 : r.l2);
        }
    }

    void push_add(int k, S x) {
        node& p = d[k];
        p.sum += x << (log + __builtin_clz(k) - 31);
        p.g1 += x;
        p.l1 += x;
        if (p.g2 != -INF) p.g2 += x;
        if (p.l2 != INF) p.l2 += x;
        p.add += x;
    }
    void push_min(int k, S x) {
        node& p = d[k];
        p.sum += (x - p.g1) * p.gc;
        if (p.l1 == p.g1) p.l1 = x;
        if (p.l2 == p.g1) p.l2 = x;
        p.g1 = x;
    }
    void push_max(int k, S x) {
        node& p = d[k];
        p.sum += (x - p.l1) * p.lc;
        if (p.g1 == p.l1) p.g1 = x;
        if (p.g2 == p.l1) p.g2 = x;
        p.l1 = x;
    }
    void push(int k) {
        node& p = d[k];
        if (p.add) {
            push_add(k * 2 + 0, p.add);
            push_add(k * 2 + 1, p.add);
            p.add = 0;
        }
        if (p.g1 < d[k * 2 + 0].g1) push_min(k * 2 + 0, p.g1);
        if (p.l1 > d[k * 2 + 0].l1) push_max(k * 2 + 0, p.l1);
    
        if (p.g1 < d[k * 2 + 1].g1) push_min(k * 2 + 1, p.g1);
        if (p.l1 > d[k * 2 + 1].l1) push_max(k * 2 + 1, p.l1);
    }

    void subtree_chmin(int k, S x) {
        if (d[k].g1 <= x) return;
        if (d[k].g2 < x) {
            push_min(k, x);
            return;
        }
        push(k);
        subtree_chmin(k * 2 + 0, x);
        subtree_chmin(k * 2 + 1, x);
        update(k);
    }

    void subtree_chmax(int k, S x) {
        if (x <= d[k].l1) return;
        if (x < d[k].l2) {
            push_max(k, x);
            return;
        }
        push(k);
        subtree_chmax(k * 2 + 0, x);
        subtree_chmax(k * 2 + 1, x);
        update(k);
    }

    template <int cmd>
    inline void _apply(int k, S x) {
        if constexpr (cmd == 1) subtree_chmin(k, x);
        if constexpr (cmd == 2) subtree_chmax(k, x);
        if constexpr (cmd == 3) push_add(k, x);
        if constexpr (cmd == 4) subtree_chmin(k, x), subtree_chmax(k, x);
    }

    template <int cmd>
    void inner_apply(int l, int r, S x) {
        if (l == r) return;
        l += size, r += size;
        for (int i = log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
        {
            int l2 = l, r2 = r;
            while (l < r) {
                if (l & 1) _apply<cmd>(l++, x);
                if (r & 1) _apply<cmd>(--r, x);
                l >>= 1;
                r >>= 1;
            }
            l = l2;
            r = r2;
        }
        for (int i = 1; i <= log; i++) {
            if (((l >> i) << i) != l) update(l >> i);
            if (((r >> i) << i) != r) update((r - 1) >> i);
        }
    }

    template <int cmd>
    inline S e() {
        if constexpr (cmd == 1) return INF;
        if constexpr (cmd == 2) return -INF;
        return 0;
    }

    template <int cmd>
    inline void op(S& a, const node& b) {
        if constexpr (cmd == 1) a = min(a, b.l1);
        if constexpr (cmd == 2) a = max(a, b.g1);
        if constexpr (cmd == 3) a += b.sum;
    }

    template <int cmd>
    S inner_fold(int l, int r) {
        if (l == r) return e<cmd>();
        l += size, r += size;
        for (int i = log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
        S lx = e<cmd>(), rx = e<cmd>();
        while (l < r) {
            if (l & 1) op<cmd>(lx, d[l++]);
            if (r & 1) op<cmd>(rx, d[--r]);
            l >>= 1;
            r >>= 1;
        }
        if constexpr (cmd == 1) lx = min(lx, rx);
        if constexpr (cmd == 2) lx = max(lx, rx);
        if constexpr (cmd == 3) lx += rx;
        return lx;
    }
};