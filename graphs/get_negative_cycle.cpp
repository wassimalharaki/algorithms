// O(VE)
vector<int> get_neg_cycle(int n, const vector<array<int, 3>>& edges) {
    vector<int> d(n), p(n, -1);
    int x;
    for (int i = 0; i < n; i++) {
        x = -1;
        for (const auto& [a, b, c] : edges)
            if (d[a] + c < d[b]) {
                x = b;
                d[x] = d[a] + c;
                p[x] = a;
            }
    }
    if (x == -1) return {};

    for (int i = 0; i < n; i++, x = p[x]);
    vector<int> cycle{x};
    for (int u = p[x]; u != x; u = p[u])
        cycle.push_back(u);
    cycle.push_back(x);
    reverse(cycle.begin(), cycle.end());

    return cycle;
}