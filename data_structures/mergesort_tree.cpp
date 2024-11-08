// O(nlog(n)), O(log^2(n))
template <class S, int (*op)(vector<S>&, S)>
struct mergesort_tree {
    int n, size, log;
    vector<vector<S>> d;

    mergesort_tree(const vector<S>& a) {
        n = a.size();
        size = n <= 1 ? 1 : 1 << (1 + __lg(n - 1));
        log = __builtin_ctz(size);
        d.resize(size << 1);
        for (int i = 0; i < n; i++) d[size + i] = {a[i]};
        for (int i = size - 1; i >= 1; i--)
            merge(
                d[i << 1].begin(), d[i << 1].end(),
                d[(i << 1) + 1].begin(), d[(i << 1) + 1].end(),
                back_inserter(d[i])
            );
    }

    S prod(int l, int r, S x) {
        l += size, r += size;

        int ans = 0;
        while (l < r) {
            if (l & 1) ans += op(d[l++], x);
            if (r & 1) ans += op(d[--r], x);
            l >>= 1;
            r >>= 1;
        }
        return ans;
    }
};