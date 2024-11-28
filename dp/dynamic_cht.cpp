struct line {
    mutable int a, b, p;
    bool operator<(const line& o) const { return a < o.a; }
    bool operator<(int x) const { return p < x; }
    int operator()(int x) const { return a * x + b; }
};

// O(log(n))
// dp[i] = max(a[j] * x[i] + b[j])
struct dynamic_cht : multiset<line, less<>> {
    static const int inf = LLONG_MAX;

    int div(int a, int b) {
        return a / b - ((a ^ b) < 0 and a % b); 
    }

    bool isect(iterator x, iterator y) {
        if (y == end())
            return x -> p = inf, 0;

        if (x -> a == y -> a)
            x -> p = x -> b > y -> b ? inf : -inf;
        else
            x -> p = div(y -> b - x -> b, x -> a - y -> a);
        return x -> p >= y -> p;
    }

    void add(int a, int b) {
        auto z = insert({a, b, 0}), y = z++, x = y;

        while (isect(y, z)) z = erase(z);
        if (x != begin() and isect(--x, y))
            isect(x, y = erase(y));

        while ((y = x) != begin() and (--x) -> p >= y -> p)
            isect(x, erase(y));
    }

    int prod(int x) { return (*lower_bound(x))(x); }
};