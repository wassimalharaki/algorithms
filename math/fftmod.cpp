using cd = complex<double>;
void fft(vector<cd>& a) {
    int n = size(a), L = 31 - __builtin_clz(n);

    static vector<complex<long double>> R(2, 1);
    static vector<cd> rt(2, 1);
    for (int k = 2; k < n; k <<= 1) {
        R.resize(n); rt.resize(n);
        auto x = polar(1.0L, acos(-1.0L) / k);
        for (int i = k; i < 2 * k; i++)
            rt[i] = R[i] = i & 1 ? R[i >> 1] * x : R[i >> 1];
    }

    vector<int> rev(n);
    for (int i = 0; i < n; i++)
        rev[i] = (rev[i >> 1] | (i & 1) << L) >> 1;

    for (int i = 0; i < n; i++)
        if (i < rev[i])
            swap(a[i], a[rev[i]]);

    for (int k = 1; k < n; k <<= 1)
        for (int i = 0; i < n; i += k << 1)
            for (int j = 0; j < k; j++) {
                auto x = (double*)& rt[j + k],
                     y = (double*)& a[i + j + k];
                cd z(x[0] * y[0] - x[1] * y[1],
                     x[0] * y[1] + x[1] * y[0]);
                a[i + j + k] = a[i + j] - z;
                a[i + j] += z;
            }
}

vector<int> mutliply(const vector<int> &a, const vector<int> &b) {
    if (a.empty() or b.empty()) return {};

    vector<int> res(size(a) + size(b) - 1);
    int B = 32 - __builtin_clz(size(res)),
        n = 1 << B, cut = (int) (sqrt(mod));

    vector<cd> L(n), R(n), outs(n), outl(n);
    for (int i = 0; i < (int) size(a); i++)
        L[i] = cd((int) a[i] / cut, (int) a[i] % cut);
    for (int i = 0; i < (int) size(b); i++)
        R[i] = cd((int) b[i] / cut, (int) b[i] % cut);
    fft(L), fft(R);

    for (int i = 0; i < n; i++) {
        int j = -i & (n - 1);
        outl[j] = (L[i] + conj(L[j])) * R[i] / (2.0 * n);
        outs[j] = (L[i] - conj(L[j])) * R[i] / (2.0 * n) / 1i;
    }
    fft(outl), fft(outs);

    for (int i = 0; i < (int) size(res); i++) {
        int av = (int) (real(outl[i]) + .5), cv = (int) (imag(outs[i]) + .5);
        int bv = (int) (imag(outl[i]) + .5) + (int) (real(outs[i]) + .5);
        res[i] = ((av % mod * cut + bv) % mod * cut + cv) % mod;
    }
    return res;
}