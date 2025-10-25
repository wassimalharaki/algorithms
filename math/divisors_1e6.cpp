const int N = 1e6 + 1;
vector<int> spf(N), primes;

// O(N)
void build() {
    for (int i = 2; i < N; i++) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (int j = 0; i * primes[j] < N; j++) {
            spf[i * primes[j]] = primes[j];
            if (primes[j] == spf[i])
                break;
        }
    }
}

// O(d(n^p))
vector<int> gen_divisors(const vector<array<int, 2>>& pfs) {
    vector<int> divs{1};

    for (auto [p, k] : pfs) {
        int n = divs.size(), x = p;
        while (k--) {
            for (int i = 0; i < n; i++)
                divs.push_back(x * divs[i]);
            x *= p;
        }
    }

    // sort(divs.begin(), divs.end());
    return divs;
}

// O(log(n) + d(n^p))
vector<int> divisors(int n, int p = 1) {
    if (primes.empty()) build();
    if (n == 1) return {1};

    vector<array<int, 2>> pfs{{spf[n], p}};
    n /= spf[n];
    while (n != 1) {
        if (pfs.back()[0] == spf[n])
            pfs.back()[1] += p;
        else
            pfs.push_back({spf[n], p});
        n /= spf[n];
    }

    return gen_divisors(pfs);
}