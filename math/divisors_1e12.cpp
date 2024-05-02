#include <bits/stdc++.h>
using namespace std;
#define int long long

// O(Nlog(log(N)))
const int N = 1e6;
bitset<N + 1> is_prime;
vector<int> primes;
void build() {
    is_prime.set();
    is_prime[0] = is_prime[1] = 0;

    primes.push_back(2);
    for (int i = 4; i <= N; i += 2) 
        is_prime[i] = 0;

    for (int i = 3; i * i <= N; i += 2)
        if (is_prime[i])
            for (int j = i * i; j <= N; j += i * 2)
                is_prime[j] = 0;

    for (int i = 3; i <= N; i += 2)
        if (is_prime[i])
            primes.push_back(i);
}

// O(d(n^p))
vector<int> gen_divisors(const vector<array<int, 2>>& pfs) {
    vector<int> divs{1};

    auto f = [&](int x, int i, auto&& f) -> void {
        if (i >= pfs.size()) return;
        f(x, i + 1, f);
        for (int j = 0; j < pfs[i][1]; j++) {
            x *= pfs[i][0];
            divs.push_back(x);
            f(x, i + 1, f);
        }
    };
    f(1, 0, f);

    // sort(divs.begin(), divs.end());
    return divs;
}

// O(sqrt(n) / log(sqrt(n)) + d(n^p))
vector<int> divisors(int n, int p = 1) {
    if (primes.empty()) build();

    vector<array<int, 2>> pfs;
    for (const int& x : primes) {
        if (x * x > n) break;
        if (n % x) continue;
        pfs.push_back({x, 0});
        while (n % x == 0) {
            n /= x;
            pfs.back()[1] += p;
        } 
    }
    if (n != 1) pfs.push_back({n, p});
    
    return gen_divisors(pfs);
}