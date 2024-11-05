// O(n^2)
// dp[i][j] = min(i <= k < j) (dp[i][k] + dp[k + 1][j] + C(i, j))
int knuth(int n, function<int(int, int)> cost) {
    vector dp(n, vector<int>(n)), opt(n, vector<int>(n));

    for (int i = 0; i < n; i++)
        opt[i][i] = i, dp[i][i] = cost(i, i);

    for (int i = n - 2; i >= 0; i--)
        for (int j = i + 1; j < n; j++) {
            int mn = INT_MAX, c = cost(i, j);
            for (int k = opt[i][j - 1]; k <= min(j - 1, opt[i + 1][j]); k++)
                if (mn >= dp[i][k] + dp[k + 1][j] + c)
                    opt[i][j] = k, mn = dp[i][k] + dp[k + 1][j] + c;
            dp[i][j] = mn;
        }
    return dp[0][n - 1];
};