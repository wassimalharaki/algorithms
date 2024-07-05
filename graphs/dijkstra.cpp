using ai2 = array<int, 2>;
// O(Elog(V))
vector<int> dijkstra(int src, vector<vector<ai2>>& adj) {
    int n = adj.size();
    vector<int> d(n, INT_MAX);
    vector<char> vis(n);
    priority_queue<ai2, vector<ai2>, greater<ai2>> pq;

    d[src] = 0;
    pq.push({0, src});

    while (pq.size()) {
        auto [w, u] = pq.top();

        if (vis[u]) continue;
        vis[u] = 1;

        for (ai2& i : adj[u])
            if (w + i[1] < d[i[0]]) {
                d[i[0]] = w + i[1];
                pq.push({d[i[0]], i[0]});
            }
    }
    return d;
}