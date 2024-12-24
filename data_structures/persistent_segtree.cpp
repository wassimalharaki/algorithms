// O(nlog(n)), O(log(n))
template <class S, S(*op)(S, S), S(*e)()>
struct pers_segtree {
    vector<int> lc, rc;
    vector<S> d;
    int n, k = 0;

    //sz >= 2(n + q) + qlog(n) should be satisfied
    pers_segtree(int sz, int _n)
        : pers_segtree(sz, vector<S>(_n, e())) {}

    pers_segtree(int sz, const vector<S>& a) {
        n = a.size();
        lc.resize(sz, -1);
        rc.resize(sz, -1);
        d.resize(sz, e());

        auto build = [&](int l, int r, auto&& self) -> int {
            int rt = k++;
            if (l + 1 == r)
                d[rt] = a[l];
            else {
                int mid = (l + r) / 2;
                lc[rt] = self(l, mid, self);
                rc[rt] = self(mid, r, self);
                d[rt] = op(d[lc[rt]], d[rc[rt]]);
            }
            return rt;
        };
        build(0, n, build);
    }

    int set(int rt, int i, S x, int l = 0, int r = -1) {
        if (r == -1) r = n;
        int nrt = k++;

        if (l + 1 == r)
            d[nrt] = x;
        else {
            if (int mid = (l + r) / 2; i < mid) {
                lc[nrt] = set(lc[rt], i, x, l, mid);
                rc[nrt] = rc[rt];
            }
            else {
                lc[nrt] = lc[rt];
                rc[nrt] = set(rc[rt], i, x, mid, r);
            }
            d[nrt] = op(d[lc[nrt]], d[rc[nrt]]);
        }
        return nrt;
    }

    S prod(int rt, int ql, int qr, int l = 0, int r = -1) {
        if (r == -1) r = n;

        if (ql == qr or rt == -1 or r <= ql or l >= qr)
            return e();
        if (ql <= l and r <= qr)
            return d[rt];
        return op(
            prod(lc[rt], ql, qr, l, (l + r) / 2),
            prod(rc[rt], ql, qr, (l + r) / 2, r)
        );
    }
};