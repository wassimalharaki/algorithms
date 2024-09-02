// O(nlog(n))
vector<int> suffix_array(string& s) {
    s.push_back('$');
    int n = s.size();
    const int alph = 256;
    vector<int> p(n), c(n), cnt(max(alph, n));

    for (int i = 0; i < n; i++)
        cnt[s[i]]++;
    for (int i = 1; i < alph; i++)
        cnt[i] += cnt[i - 1];
    for (int i = 0; i < n; i++)
        p[--cnt[s[i]]] = i;

    c[p[0]] = 0;
    int cls = 1;
    for (int i = 1; i < n; i++) {
        if (s[p[i]] != s[p[i - 1]])
            cls++;
        c[p[i]] = cls - 1;
    }

    vector<int> pn(n), cn(n);
    for (int h = 0; (1 << h) < n; h++) {
        for (int i = 0; i < n; i++) {
            pn[i] = p[i] - (1 << h);
            if (pn[i] < 0) pn[i] += n;
        }
        fill(cnt.begin(), cnt.begin() + cls, 0);

        for (int i = 0; i < n; i++)
            cnt[c[pn[i]]]++;
        for (int i = 1; i < cls; i++)
            cnt[i] += cnt[i - 1];
        for (int i = n - 1; i >= 0; i--)
            p[--cnt[c[pn[i]]]] = pn[i];

        cn[p[0]] = 0;
        cls = 1;
        for (int i = 1; i < n; i++) {
            array<int, 2> curr = {
                c[p[i]],
                c[(p[i] + (1 << h)) % n]
            };
            array<int, 2> prev = {
                c[p[i - 1]],
                c[(p[i - 1] + (1 << h)) % n]
            };
            if (curr != prev) cls++;
            cn[p[i]] = cls - 1;
        }
        c.swap(cn);
    }
    
    p.erase(p.begin());
    s.pop_back();
    return p;
}