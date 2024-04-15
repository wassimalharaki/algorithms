#include <bits/stdc++.h>
using namespace std;

// O(N)
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