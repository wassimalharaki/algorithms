// O(n)
int f(string& s, int i = 0, bool t = 1) {
    if (dp[i][t] != -1)
        return dp[i][t];
    int& ans = dp[i][t] = 0;
    int r = t ? s[i] - '0' : 9;

    for (int j = 0; j <= r; j++)
        ans += f(s, i + 1, x - 1, t and i == r);
    return ans;
}