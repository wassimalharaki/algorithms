// O(n / 64), O(log_64(n))
struct fast_set {
    using u64 = uint64_t;
    static constexpr uint32_t B = 64;
    int n, log;
    vector<vector<u64>> seg;

    fast_set() {}
    fast_set(int _n) { build(_n); }

    void build(int m) {
        seg.clear();
        n = m;
        do {
            seg.push_back(vector<u64>((m + B - 1) / B));
            m = (m + B - 1) / B;
        } while (m > 1);
        log = seg.size();
    }

    bool operator[](int i) const {
        return seg[0][i / B] >> (i % B) & 1;
    }
    void insert(int i) {
        for (int h = 0; h < log; h++)
            seg[h][i / B] |= u64(1) << (i % B), i /= B;
    }
    void erase(int i) {
        u64 x = 0;
        for (int h = 0; h < log; h++) {
            seg[h][i / B] &= ~(u64(1) << (i % B));
            seg[h][i / B] |= x << (i % B);
            x = bool(seg[h][i / B]);
            i /= B;
        }
    }

    // min[x, n) or n
    int next(int i) {
        i = max(i, 0);
        for (int h = 0; h < log; h++) {
            if (i / B == seg[h].size()) break;
            u64 d = seg[h][i / B] >> (i % B);
            if (!d) {
                i = i / B + 1;
                continue;
            }
            i += __builtin_ctzll(d);
            for (int g = h - 1; g >= 0; g--) {
                i *= B;
                i += __builtin_ctzll(seg[g][i / B]);
            }
            return i;
        }
        return n;
    }

    // max[0, x] or -1
    int prev(int i) {
        if (i >= n) i = n - 1;
        for (int h = 0; h < log; h++) {
            if (i == -1) break;
            u64 d = seg[h][i / B] << (63 - i % B);
            if (!d) {
                i = i / B - 1;
                continue;
            }
            i -= __builtin_clzll(d);
            for (int g = h - 1; g >= 0; g--) {
                i *= B;
                i += __lg(seg[g][i / B]);
            }
            return i;
        }
        return -1;
    }
};