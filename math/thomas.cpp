// O(nlog(b[i]))
// a[i] * f(x - 1) + b[i] * f(x) + c[i] * f(x + 1) = d[i]
vector<int> thomas(vector<int>& a, vector<int>& b, vector<int>& c, vector<int>& d) {
    int n = a.size();
    vector<int> x(n);
    for (int i = 1; i < n; i++) {
        int w = a[i] * modinv(b[i - 1]) % mod;
        b[i] = (b[i] - w * c[i - 1] % mod + mod) % mod;
        d[i] = (d[i] - w * d[i - 1] % mod + mod) % mod;
    }
    x[n - 1] = d[n - 1] * modinv(b[n - 1]) % mod;
    for (int i = n - 2; i >= 0; i--)
        x[i] = (d[i] - c[i] * x[i + 1] % mod + mod) * modinv(b[i]) % mod;
    return x;
}