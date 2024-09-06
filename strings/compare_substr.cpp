// O(nlog(n)), O(1)
struct compare_substr {
    vector<vector<int>> c;

    compare_substr(string& s) {
        s.push_back('$');
        int n = s.size(), k = 2 + (n == 1 ? 0 : __lg(n - 1));
        const int alph = 256;
        vector<int> p(n), cnt(max(alph, n));
        c.resize(k, vector<int>(n));

        for (int i = 0; i < n; i++)
            cnt[s[i]]++;
        for (int i = 1; i < alph; i++)
            cnt[i] += cnt[i - 1];
        for (int i = 0; i < n; i++)
            p[--cnt[s[i]]] = i;

        c[0][p[0]] = 0;
        int cls = 1;
        for (int i = 1; i < n; i++) {
            if (s[p[i]] != s[p[i - 1]])
                cls++;
            c[0][p[i]] = cls - 1;
        }

        vector<int> pn(n);
        for (int j = 1; j < k; j++) {
            for (int i = 0; i < n; i++) {
                pn[i] = p[i] - (1 << (j - 1));
                if (pn[i] < 0) pn[i] += n;
            }
            fill(cnt.begin(), cnt.begin() + cls, 0);

            for (int i = 0; i < n; i++)
                cnt[c[j - 1][pn[i]]]++;
            for (int i = 1; i < cls; i++)
                cnt[i] += cnt[i - 1];
            for (int i = n - 1; i >= 0; i--)
                p[--cnt[c[j - 1][pn[i]]]] = pn[i];

            c[j][p[0]] = 0, cls = 1;
            for (int i = 1; i < n; i++) {
                array<int, 2> cur = {
                    c[j - 1][p[i]],
                    c[j - 1][(p[i] + (1 << (j - 1))) % n]
                };
                array<int, 2> prev = {
                    c[j - 1][p[i - 1]],
                    c[j - 1][(p[i - 1] + (1 << (j - 1))) % n]
                };
                if (cur != prev) cls++;
                c[j][p[i]] = cls - 1;
            }
        }
        s.pop_back();
    }

    int prod(int l1, int r1, int l2, int r2) {
        if (r1 - l1 > r2 - l2)
            return - prod(l2, r2, l1, r1);

        if (l1 == r1)
            return l2 == r2 ? 0 : -1;

        if (r1 - l1 != r2 - l2) {
            int x = prod(l1, r1, l2, l2 + r1 - l1);
            return x ? x : -1;
        }

        int i = __lg(r1 - l1);
        array<int, 2> a = {c[i][l1], c[i][r1 - (1 << i)]};
        array<int, 2> b = {c[i][l2], c[i][r2 - (1 << i)]};
        return a == b ? 0 : a < b ? -1 : 1;
    }
};