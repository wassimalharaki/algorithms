template <class S, class F>
struct mo_s {
    struct query { int l, r, i; int64_t h; };
    vector<query> q;
    int m = 0;

    int64_t hilbert_order(int l, int r){
        int lg = __lg((r << 1) + 1) | 1;
        int mxN = (1 << lg) - 1;
        int64_t ans = 0;
        for (int i = 1 << (lg - 1); i; i >>= 1) {
            bool rl = l & i, ry = r & i;
            ans = (ans << 2) | (rl ? ry ? 2 : 1 : ry ? 3 : 0);
            if (not rl) {
                if (ry) l ^= mxN, r ^= mxN;
                swap(l, r);
            }
        }
        return ans;
    }

    void prod(int l, int r) {
        q.push_back({l, r, m++, hilbert_order(l, r)});
    }

    vector<F> solve(const vector<S>& a) {
        sort(q.begin(), q.end(), [&](auto& x, auto& y) {
            return x.h < y.h;
        });

        auto insert = [&](int i, bool back) {

        };

        auto erase = [&](int i, bool back) {

        };

        int l = 0, r = 0;
        vector<F> ans(m);
        for (auto& [_l, _r, i, _] : q) {
            while (r < _r) insert(r++, 1);
            while (l > _l) insert(--l, 0);
            while (r > _r) erase(--r, 1);
            while (l < _l) erase(l++, 0);
            ans[i] = 0;
        }
        return ans;
    }
};