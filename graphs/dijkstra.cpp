#include <bits/stdc++.h>
using namespace std;

#define v vector
#define F first
#define S second

using pii = pair<int, int>;

// O(Elog(V))
v<int> dijkstra(int src, v<v<pii>>& adj) {
    int n = adj.size();
    v<int> d(n, INT_MAX);
    v<bool> vis(n);
    priority_queue<pii, v<pii>, greater<pii>> pq;

    d[src] = 0;
    pq.push({0, src});

    while (pq.size()) {
        auto [w, u] = pq.top();

        if (vis[u]) continue;
        vis[u] = 1;

        for (pii& i : adj[u])
            if (w + i.F < d[i.S]) {
                d[i.S] = w + i.F;
                pq.push({d[i.S], i.S});
            }
    }
    return d;
}