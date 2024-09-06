const int N = 5e5 + 1;
const int mod = 1e9 + 7;
int part[N];

int pent(int n) {
    return n * (3 * n - 1) >> 1;
}

void add(int& a, int b) {
    a += b;
    while (a < 0) a += mod;
    while (a >= mod) a -= mod;
}

// O(Nsqrt(N))
void build() {
    part[0] = 1;
    for (int i = 1; i < N; i++)
        for (int j = 1; j <= i and i - pent(j) >= 0; j++)
            for (int k : {pent(j), pent(-j)})
                if (i - k >= 0)
                    add(part[i], (j & 1 ? 1 : -1) * part[i - k]);
}