int dp[20][2];

int f(string& s, int i = 0, bool t = 1) {
    if (dp[n][t] != -1)
        return dp[n][t];
    int& ans = dp[n][t] = 0;
    int r = t ? s[i] - '0' : 9;

    for (int j = 0; j <= r; j++)
        ans += f(s, i + 1, x - 1, t and i == r);
    return ans;
}