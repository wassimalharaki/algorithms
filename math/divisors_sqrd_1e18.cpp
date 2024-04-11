#include <bits/stdc++.h>
using namespace std;
#define int long long

const int N = 4e4;
bitset<N + 1> is_prime;
vector<int> primes{2};
void build() {
    is_prime.set();
    is_prime[0] = is_prime[1] = 0;

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

// O(sqrt(n) / log(sqrt(n)) + d(n^2))
vector<int> divisors_sqrd(int n) {
    if (primes.size() == 1) build();

    vector<pair<int, int>> pfs;
    for (int& p : primes) {
        if (p * p > n) break;
        if (n % p) continue;
        pfs.push_back({p, 0});
        while (n % p == 0) {
            n /= p;
            pfs.back().second += 2;
        } 
    }
    if (n != 1) pfs.push_back({n, 2});
    
    return gen_divisors(pfs);
}