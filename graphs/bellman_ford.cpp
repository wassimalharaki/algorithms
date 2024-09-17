// O(VE)
vector<int> bellman_ford(int src, int n, const vector<vector<int>>& edges) {
    vector<int> d(n, INT_MAX);
    d[src] = 0;

    for (int i = 0; i < n - 1; i++)
        for (const auto& [a, b, c] : edges)
            if (d[a] + c < d[b])
                d[b] = d[a] + c;

    for (const auto& [a, b, c] : edges)
        if (d[a] + c < d[b])
            return {};

    return d;
}