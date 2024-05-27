#include <bits/stdc++.h>
using namespace std;

// O(Nlog(log(N)))
const int N = 1e6 + 1;
bitset<N> is_prime;
vector<int> primes{2};
void build() {
    is_prime.set();
    is_prime[0] = is_prime[1] = 0;

    for (int i = 4; i < N; i += 2) 
        is_prime[i] = 0;

    for (int i = 3; i * i < N; i += 2)
        if (is_prime[i])
            for (int j = i * i; j < N; j += i * 2)
                is_prime[j] = 0;

    for (int i = 3; i < N; i += 2)
        if (is_prime[i])
            primes.push_back(i);
}