// O(n), O(log(n))
template <class T>
struct BIT {
    int n;
    vector<T> d;

    BIT(int _n) {
        n = _n;
        d.resize(_n);
    }

    T sum(int r) {
        T s = 0;
        while (r > 0) {
            s += d[r - 1];
            r -= r & -r;
        }
        return s;
    }

    void add(int p, T x) {
        p++;
        while (p <= n) {
            d[p - 1] += x;
            p += p & -p;
        }
    }

    T sum(int l, int r) {
        return sum(r) - sum(l);
    }
};