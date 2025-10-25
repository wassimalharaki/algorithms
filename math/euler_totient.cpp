const int N = 1e6 + 1;
vector<int> phi(N);

// O(Nlog(log(N)))
void build() {
    iota(phi.begin(), phi.end(), 0);
    for (int i = 2; i < N; i++)
        if (phi[i] == i)
            for (int j = i; j < N; j += i)
                phi[j] -= phi[j] / i;
}