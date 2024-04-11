#include <bits/stdc++.h>
using namespace std;

struct segtree {
    int n, size, log;
    vector<int> d;

    segtree(int n) : segtree(vector<int>(n, 0)) {}

    segtree(const vector<int>& nums) {
        n = nums.size();
        size = bit_ceil(n);
        log = __builtin_ctz(size);
        d = vector<int>(2 * size, 0);
        for (int i = 0; i < n; i++) d[size + i] = nums[i];
        for (int i = size - 1; i >= 1; i--) update(i);
    }

    void update(int k) { d[k] = d[2 * k] + d[2 * k + 1]; }
    
    int bit_ceil(int n) {
        int x = 1;
        while (x < n) x <<= 1;
        return x;
    }

    void set(int p, int x) {
        p += size;
        d[p] = x;
        for (int i = 1; i <= log; i++) update(p >> i);
    }

    int get(int p) const {
        return d[p + size];
    }

    int sum(int l, int r) const {
        int sml = 0, smr = 0;
        l += size;
        r += size;

        while (l < r) {
            if (l & 1) sml = sml + d[l++];
            if (r & 1) smr = d[--r] + smr;
            l >>= 1;
            r >>= 1;
        }
        return sml + smr;
    }

    int all_sum() const { return d[1]; }
};