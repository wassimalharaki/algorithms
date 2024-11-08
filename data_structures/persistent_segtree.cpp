template <class S, S(*op)(S, S), S(*e)()>
struct pers_segtree {
    vector<int> lc, rc;
    vector<S> d, a;
    int n = 0;

    //sz >= (2n - 1) + q * (bit_width(n) + 1) should be satisfied
    pers_segtree(int sz, vector<S> _a = vector<S>())
    : lc(sz, -1), rc(sz, -1), d(sz, e()), a(_a) {}

    int build(int l, int r) {
        int i = n++;
        if (l + 1 == r) {
            if (a.size())
                d[i] = a[l];
        }
        else {
            int mid = (l + r) / 2;
            lc[i] = build(l, mid), rc[i] = build(mid, r);
            d[i] = op(d[lc[i]], d[rc[i]]);
        }
        return i;
    }

    int set(int rt, int l, int r, int i, S x) {  
        int nrt = n++;
        if (l + 1 == r)
            d[nrt] = x;
        else {
            if (int mid = (l + r) / 2; i < mid)
                lc[nrt] = set(lc[rt], l, mid, i, x), rc[nrt] = rc[rt];
            else
                lc[nrt] = lc[rt], rc[nrt] = set(rc[rt], mid, r, i, x);
            d[nrt] = op(d[lc[nrt]], d[rc[nrt]]);
        }
        return nrt;
    }

    S prod(int i, int l, int r, int ql, int qr) {
        if (ql == qr or i == -1 or r <= ql or l >= qr)
            return e();
        if (ql <= l and r <= qr)
            return d[i];
        return op(
            prod(lc[i], l, (l + r) / 2, ql, qr), 
            prod(rc[i], (l + r) / 2, r, ql, qr)
        );
    }
};