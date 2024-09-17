// O(n * m)
template <class S>
struct prefix_2d {
    vector<vector<S>> d;

    prefix_2d(const vector<vector<S>>& a) {
        int n = a.size(), m = a[0].size();
        d.resize(n + 1, vector<S>(m + 1));

        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= m; j++)
                d[i][j] = d[i - 1][j] + d[i][j - 1]
                    + a[i - 1][j - 1] - d[i - 1][j - 1];
    }

    S prod(int i1, int j1, int i2, int j2) {
        return d[i2][j2] + d[i1][j1]
            - d[i2][j1] - d[i1][j2];
    }
};