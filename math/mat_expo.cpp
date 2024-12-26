// O(nm), O(nm^2)
template <class S, S mod = (S) 1e9 + 7>
struct matrix {
    vector<vector<S>> a;

    matrix(const vector<vector<S>>& _a) { a = _a; }
    matrix(int n, int m) { a.resize(n, vector<S>(m)); }

    size_t size() const { return a.size(); }
    vector<S>& operator[](int i) { return a[i]; }
    const vector<S>& operator[](int i) const { return a[i]; }

    matrix& operator*=(const matrix& b) {
        int ra = a.size(), ca = a[0].size(), cb = b[0].size();
        matrix c(ra, cb);
        for (int i = 0; i < ra; i++)
            for (int j = 0; j < cb; j++)
                for (int k = 0; k < ca; k++)
                    c[i][j] = (c[i][j] + a[i][k] * b[k][j]) % mod;
        a.swap(c.a);
        return *this;
    }
};

// O(n^3 log(b))
matrix<int> binpow(matrix<int> a, int b) {
    int n = a.size();
    matrix<int> res(n, n);
    for (int i = 0; i < n; i++)
        res[i][i] = 1;

    while (b) {
        if (b & 1) res *= a;
        a *= a, b >>= 1;
    }
    return res;
}