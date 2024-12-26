// O(N), O(N^2)
template <class S, int N>
struct basis {
    array<S, N> a;
    int n;

    basis() { a.fill(0); n = 0; }
    void clear() { a.fill(0); n = 0; }

    void insert(S x) {
        for (int i = N - 1; i >= 0; i--)
            x = min(x, x ^ a[i]);
        if (x) a[__lg(x)] = x, n++;
    }

    void merge(const basis<S, N>& b) {
        for (int x : b.a) insert(x);
    }
};