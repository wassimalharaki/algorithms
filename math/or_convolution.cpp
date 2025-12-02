// O(nlog(n))
vector<int> or_convolution(vector<int> a, vector<int> b) {
    int m = __lg(a.size());

    for (int i = 0; i < m; i++)
        for (int j = 0; j < (1 << m); j++)
            if ((1 << i) & j) {
                a[j] += a[j ^ (1 << i)];
                b[j] += b[j ^ (1 << i)];
            }

    for (int i = 0; i < (1 << m); i++)
        a[i] *= b[i];

    for (int i = 0; i < m; i++)
        for (int j = 0; j < (1 << m); j++)
            if ((1 << i) & j)
                a[j] -= a[j ^ (1 << i)];
    return a;
}