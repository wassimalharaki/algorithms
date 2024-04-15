#include <bits/stdc++.h>
using namespace std;

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

// O(sqrt(n) / log(sqrt(n)))
vector<array<int, 2>> prime_factors(int n) {
    if (primes.empty()) build();

    vector<array<int, 2>> pfs;
    for (const int& x : primes) {
        if (x * x > n) break;
        if (n % x) continue;
        pfs.push_back({x, 0});
        while (n % x == 0) {
            n /= x;
            pfs.back()[1]++;
        } 
    }
    if (n != 1) pfs.push_back({n, 1});
    
    return pfs;
}