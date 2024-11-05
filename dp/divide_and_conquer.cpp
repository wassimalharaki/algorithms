// O(nmlog(n))
// dp[i][j] = min(0 <= k <= j) (dp[i - 1][k - 1] + C(k, j))
vector<int> dnc(int n, int m, function<int(int, int)> cost) {
    vector<int> dp(n), ndp(n);
    auto f = [&](int l, int r, int optl, int optr, auto&& self) -> void {
        if (l > r) return;
        int mid = (l + r) >> 1;
        array<int, 2> best{LLONG_MAX, -1};
        
        for (int k = optl; k <= min(mid, optr); k++)
            best = min(best, {(k ? dp[k - 1] : 0) + cost(k, mid), k});

        ndp[mid] = best[0];
        int opt = best[1];
        
        self(l, mid - 1, optl, opt, self);
        self(mid + 1, r, opt, optr, self);
    };

    for (int i = 0; i < n; i++)
        dp[i] = cost(0, i);
    
    for (int i = 1; i < m; i++) {
        f(0, n - 1, 0, n - 1, f);
        swap(dp, ndp);
    }
    return dp;
}