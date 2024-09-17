using ai2 = array<int, 2>;
// O(Elog(V))
vector<int> dijkstra(int src, const vector<vector<ai2>>& adj) {
    int n = adj.size();
    vector<int> d(n, INT_MAX);
    vector<char> vis(n);
    priority_queue<ai2, vector<ai2>, greater<ai2>> pq;

    d[src] = 0;
    pq.push({0, src});

    while (pq.size()) {
        auto [w, u] = pq.top();
        pq.pop();

        if (vis[u]) continue;
        vis[u] = 1;

        for (const auto& [i, c] : adj[u])
            if (w + c < d[i]) {
                d[i] = w + c;
                pq.push({d[i], i});
            }
    }
    return d;
}