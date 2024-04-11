#include <bits/stdc++.h>
using namespace std;

struct lazy_segtree {
    int n, size, log;
    vector<int> d, lz;

    lazy_segtree(int n) : lazy_segtree(vector<int>(n)) {}

    lazy_segtree(const vector<int>& nums) {
        n = nums.size();
        size = bit_ceil(n);
        log = __builtin_ctzll(size);
        d = vector<int>(size << 1);
        lz = vector<int>(size);
        for (int i = 0; i < n; i++) d[size + i] = nums[i];
        for (int i = size - 1; i >= 1; i--) update(i);
    }

    void update(int k) { d[k] = d[k << 1] + d[(k << 1) + 1]; }

    void all_apply(int k, int f) {
        d[k] = f + d[k];
        if (k < size) lz[k] = f + lz[k];
    }

    void push(int k) {
        all_apply(k << 1, lz[k]);
        all_apply((k << 1) + 1, lz[k]);
        lz[k] = 0;
    }

    int bit_ceil(int n) {
        int x = 1;
        while (x < n) x <<= 1;
        return x;
    }

    void set(int p, int x) {
        p += size;
        for (int i = log; i >= 1; i--) push(p >> i);
        d[p] = x;
        for (int i = 1; i <= log; i++) update(p >> i);
    }

    int get(int p) {
        p += size;
        for (int i = log; i >= 1; i--) push(p >> i);
        return d[p];
    }

    int sum(int l, int r) {
        if (l == r) return 0;

        l += size;
        r += size;

        for (int i = log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }

        int sml = 0, smr = 0;
        while (l < r) {
            if (l & 1) sml = sml + d[l++];
            if (r & 1) smr = d[--r] + smr;
            l >>= 1;
            r >>= 1;
        }

        return sml + smr;
    }

    int all_sum() { return d[1]; }

    void apply(int p, int f) {
        p += size;
        for (int i = log; i >= 1; i--) push(p >> i);
        d[p] = f + d[p];
        for (int i = 1; i <= log; i++) update(p >> i);
    }

    void apply(int l, int r, int f) {
        if (l == r) return;

        l += size;
        r += size;

        for (int i = log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }

        {
            int l2 = l, r2 = r;
            while (l < r) {
                if (l & 1) all_apply(l++, f);
                if (r & 1) all_apply(--r, f);
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
};