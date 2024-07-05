// O(N)
const int N = 1e6 + 1;
vector<int> spf(N), mobius(N), primes;
void build() {
    mobius[1] = 1;
    for (int i = 2; i < N; i++) {
        if (spf[i] == 0) {
            spf[i] = i;
            mobius[i] = -1;
            primes.push_back(i);
        }
        for (int j = 0; i * primes[j] < N; j++) {
            spf[i * primes[j]] = primes[j];
            if (i % primes[j])
                mobius[i * primes[j]] = mobius[i] * mobius[primes[j]];
            else
                mobius[i * primes[j]] = 0;
            if (primes[j] == spf[i])
                break;
        }
    }
}