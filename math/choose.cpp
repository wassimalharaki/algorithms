const int mod = 1e9 + 7;
const int N = 1e6 + 1;
int fact[N], inv_num[N], inv_fact[N];

void build() {
    fact[0] = fact[1] = 1;
    inv_num[0] = inv_num[1] = 1;
    inv_fact[0] = inv_fact[1] = 1;
    for (int i = 2; i < N; i++) {
        fact[i] = i * fact[i - 1] % mod;
        inv_num[i] = (mod - mod / i) * inv_num[mod % i] % mod;
        inv_fact[i] = inv_fact[i - 1] * inv_num[i] % mod;
    }
}

int choose(int n, int r) {
    return fact[n] * inv_fact[r] % mod * inv_fact[n - r] % mod;
}