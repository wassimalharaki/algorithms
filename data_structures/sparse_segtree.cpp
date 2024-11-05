// O(log(n))
template <class F, class S, S (*op)(S, S), S (*e)()>
struct sparse_segtree {
    struct node {
        F l, r;
        S x = e();
        int lc = 0, rc = 0, p;
    };
    vector<node> d;

    int add_node(F l, F r, int p) {
        d.push_back(node());
        tie(d.back().l, d.back().r, d.back().p) = {l, r, p};
        return d.size() - 1;
    }

    int get_lc(int i) {
        if (!d[i].lc)
            d[i].lc = add_node(d[i].l, (d[i].l + d[i].r) >> 1, i);
        return d[i].lc;
    }

    int get_rc(int i) {
        if (!d[i].rc)
            d[i].rc = add_node((d[i].l + d[i].r) >> 1, d[i].r, i);
        return d[i].rc;
    }

    void update(int i) {
        d[i].x = op(d[d[i].lc].x, d[d[i].rc].x);
    }

    sparse_segtree(F l, F r) {
        add_node((F) -1, (F) -1, -1);
        add_node(l, r, 0);
    }

    S get(F i) {
        int j = 1;
        while (j and (i != d[j].l or d[j].l + 1 != d[j].r))
            j = i < (d[j].l + d[j].r) >> 1 ?
                d[j].lc :
                d[j].rc;
        return d[j].x;
    }

    void set(F i, S x) {
        int j = 1;
        while (d[j].l != i or d[j].l + 1 != d[j].r)
            j = i < (d[j].l + d[j].r) >> 1 ?
                get_lc(j) :
                get_rc(j);

        d[j].x = x;
        while (d[j].p) {
            j = d[j].p;
            update(j);
        }
    }

    void apply(F i, S x) {
        int j = 1;
        d[j].x = op(d[j].x, x);
        while (d[j].l != i or d[j].l + 1 != d[j].r) {
            j = i < (d[j].l + d[j].r) >> 1 ?
                get_lc(j) :
                get_rc(j);
            d[j].x = op(d[j].x, x);
        }
    }

    S prod(F l, F r, int j = 1) {
        if (j == 0 or d[j].r <= l or d[j].l >= r)
            return e();
        if (l <= d[j].l and d[j].r <= r)
            return d[j].x;
        return op(
            prod(l, r, d[j].lc),
            prod(l, r, d[j].rc)
        );
    }
};