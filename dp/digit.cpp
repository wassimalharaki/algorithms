int dp[20][2][2];

int f(string& s, int n, bool lz = 1, bool t = 1) {
    if (dp[n][lz][t] != -1)
        return dp[n][lz][t];
    int& ans = dp[n][lz][t] = 0;

    int ans = 0;
    for (int i = 0; i <= r; i++)
        if (i == 0)
            ans += f(s, n - 1, x, lz, 0);
        else
            ans += f(s, n - 1, x - 1, 0, 0);

    return ans;
}