#include <bits/stdc++.h>
using namespace std;

#define v vector
#define F first
#define S second

using pii = pair<int, int>;

v<int> dijkstra(int src, v<v<pii>>& adj) {
    v<int> d(adj.size(), INT_MAX);
    priority_queue<pii, v<pii>, greater<pii>> pq;
    v<int> vis(adj.size(), 0);

    d[src] = 0;
    pq.push({0, src});

    while (not pq.empty()) {
        int u = pq.top().S;
        int w = pq.top().F;
        pq.pop();

        if (vis[u]) continue;
        vis[u] = 1;

        for (auto i : adj[u])
            if (w + i.F < d[i.S]) {
                d[i.S] = w + i.F;
                pq.push({d[i.S], i.S});
            }
    }
    return d;
}