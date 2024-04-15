#include <bits/stdc++.h>
using namespace std;

// O(nlog(n)), O(log(n))
template <class S, S (*op)(S, S), S (*e)()>
struct segtree {
    int n, size, log;
    vector<S> d;

    segtree(int n) : segtree(vector<S>(n, e())) {}
    
    segtree(const vector<S>& nums) {
        n = nums.size();
        size = bit_ceil(n);
        log = __builtin_ctz(size);
        d = vector<S>(size << 1, e());
        for (int i = 0; i < n; i++) d[size + i] = nums[i];
        for (int i = size - 1; i >= 1; i--) update(i);
    }

    void update(int k) { d[k] = op(d[k << 1], d[(k << 1) + 1]); }

    int bit_ceil(int n) {
        int x = 1;
        while (x < n) x <<= 1;
        return x;
    }

    void set(int p, S x) {
        p += size;
        d[p] = x;
        for (int i = 1; i <= log; i++) update(p >> i);
    }

    S get(int p) const {
        return d[p + size];
    }

    S prod(int l, int r) const {
        S sml = e(), smr = e();
        l += size;
        r += size;

        while (l < r) {
            if (l & 1) sml = op(sml, d[l++]);
            if (r & 1) smr = op(d[--r], smr);
            l >>= 1;
            r >>= 1;
        }
        return op(sml, smr);
    }

    S all_prod() const { return d[1]; }
};