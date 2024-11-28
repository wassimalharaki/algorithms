#include "./polygon.cpp"

ftype mindist, minper, mxdist;
pair<P, P> bp, bp_mx;
vec<P> btr;
int n;
vec<P> a, t;

bool cmp_x(const P& p0, const P& p1) {
    return p0 < p1;
}

bool cmp_y(const P& p0, const P& p1) {
    return ls(p0.y, p1.y);
}

void upd_ans_t(const P& p0, const P& p1, const P& p2) {
    ftype per = p0.dist(p1) + p0.dist(p2) + p1.dist(p2);
    if (ls(per, minper)) {
        minper = per;
        btr = {p0, p1, p2};
    }
}

void upd_ans(const P& p0, const P& p1) {
    ftype d = p0.dist(p1);
    if (ls(d, mindist)) {
        mindist = d;
        bp = {p0, p1};
    }
}

void rec_t(int l, int r) {
    if (r - l <= 3) {
        // Check all triangles in this small range
        for (int i = l; i < r; ++i) {
            for (int j = i + 1; j < r; ++j) {
                for (int k = j + 1; k < r; ++k) {
                    upd_ans_t(a[i], a[j], a[k]);
                }
            }
        }
        sort(a.begin() + l, a.begin() + r, cmp_y);
        return;
    }

    int m = (l + r) >> 1;
    ftype midx = a[m].x;
    rec_t(l, m);
    rec_t(m, r);

    merge(a.begin() + l, a.begin() + m, a.begin() + m, a.begin() + r, t.begin(), cmp_y);
    copy(t.begin(), t.begin() + r - l, a.begin() + l);

    vec<P> strip;
    for (int i = l; i < r; ++i) {
        if (ls(abs(a[i].x - midx), minper / 2)) {
            for (int j = (int)strip.size() - 1; j >= 0 && ls(a[i].y - strip[j].y, minper / 2); --j) {
                for (int k = j - 1; k >= 0 && ls(a[i].y - strip[k].y, minper / 2); --k) {
                    upd_ans_t(a[i], strip[j], strip[k]);
                }
            }
            strip.push_back(a[i]);
        }
    }
}

void rec(int l, int r) {
    if (r - l <= 3) {
        for (int i = l; i < r; ++i) {
            for (int j = i + 1; j < r; ++j) {
                upd_ans(a[i], a[j]);
            }
        }
        sort(a.begin() + l, a.begin() + r, cmp_y);
        return;
    }

    int m = (l + r) >> 1;
    ftype midx = a[m].x;
    rec(l, m);
    rec(m, r);

    merge(a.begin() + l, a.begin() + m, a.begin() + m, a.begin() + r, t.begin(), cmp_y);
    copy(t.begin(), t.begin() + r - l, a.begin() + l);

    int tsz = 0;
    for (int i = l; i < r; ++i) {
        if (ls(abs(a[i].x - midx), mindist)) {
            for (int j = tsz - 1; j >= 0 and ls(a[i].y - t[j].y, mindist); --j) upd_ans(a[i], t[j]);
            t[tsz++] = a[i];
        }
    }
}

void find_mn() {
    t.resize(n);
    sort(all(a), cmp_x);
    mindist = 1e20;
    rec(0, n);
}

void find_mn_t() {
    t.resize(n);
    sort(all(a), cmp_x);
    minper = 1e20;
    rec_t(0, n);
}

void find_mx() {
    // rotating caliper
    mxdist = 0;
    auto ch = convexHullGrahms(a).ps;
    reverse(all(ch)); // change from ccw to cw
    int m = (int)ch.size();
    int k = 1;
    while (k < m) {
        ftype cur = heron(ch[0], ch[m - 1], ch[k]);
        ftype nxt = heron(ch[0], ch[m - 1], ch[k + 1]);
        if (ls(nxt, cur)) break;
        k++;
    }
    int p = 0, q = k;
    while (p <= k and q < m) {
        while (q < m) {
            ftype cur = heron(ch[p], ch[p + 1], ch[q]);
            ftype nxt = heron(ch[p], ch[p + 1], ch[(q + 1) % m]);
            if (ls(nxt, cur)) break;
            q++;
            ftype tmp = ch[p].dist(ch[q % m]);
            if (ls(mxdist, tmp)) {
                mxdist = tmp;
                bp_mx = {ch[p], ch[q % m]};
            }
        }
        p++;
    }
}
