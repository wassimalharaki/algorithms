#include <bits/stdc++.h>
using namespace std;
#define int long long

const int N = 1e6;
vector<int> spf, primes;
void build() {
    spf = vector<int>(N + 1);
    for (int i = 2; i <= N; i++) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (int j = 0; i * primes[j] <= N; j++) {
            spf[i * primes[j]] = primes[j];
            if (primes[j] == spf[i])
                break;
        }
    }
}

vector<int> gen_divisors(vector<pair<int, int>>& pfs) {
    vector<int> divs{1};
    auto f = [&](int x, int i, auto&& f) -> void {
        if (i >= pfs.size()) return;
        f(x, i + 1, f);
        for (int j = 0; j < pfs[i].second; j++) {
            x *= pfs[i].first;
            divs.push_back(x);
            f(x, i + 1, f);
        }
    };
    f(1, 0, f);

    // sort(divs.begin(), divs.end());
    return divs;
}

// O(log(n) + d(n^2))
vector<int> divisors_sqrd(int n) {
    if (n == 1) return {1};
    if (primes.empty()) build();

    vector<pair<int, int>> pfs{{spf[n], 2}};
    n /= spf[n];
    while (n != 1) {
        if (pfs.back().first == spf[n])
            pfs.back().second += 2;
        else
            pfs.push_back({spf[n], 2});
        n /= spf[n];
    }

    return gen_divisors(pfs);
}