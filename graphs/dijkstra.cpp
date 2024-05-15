#include <bits/stdc++.h>
using namespace std;
#define v vector

using ai2 = array<int, 2>;

// O(Elog(V))
v<int> dijkstra(int src, v<v<ai2>>& adj) {
    int n = adj.size();
    v<int> d(n, INT_MAX);
    v<char> vis(n);
    priority_queue<ai2, v<ai2>, greater<ai2>> pq;

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