#include <bits/stdc++.h>
using namespace std;
#define v vector

// O(n), O(log(n))
struct DSU {
    vector<int> p;

    DSU(int n) { p = vector<int>(n, -1); }

    int find(int x) { return p[x] < 0 ? x : p[x] = find(p[x]); }

    int size(int x) { return - p[find(x)]; }

    bool same_set(int x, int y) { return find(x) == find(y); }

    bool merge(int x, int y) {
        x = find(x); y = find(y);
        if (x == y) return false;
        if (p[x] > p[y]) swap(x, y);
        p[x] += p[y]; p[y] = x;
        return true;
    }

    int components() {
        int x = 0;
        for (int& i : p) if (i < 0) x++;
        return x;
    }
};

// O(nlog(n))
int kruskal(int n, v<array<int, 3>>& edges) {
    sort(edges.begin(), edges.end());

    // v<array<int, 3>> tree(n - 1, array<int, 3>(3));
    DSU ds(n + 1);
    int cost = 0, j = 0;
    for (int i = 0; j < n - 1 and i < edges.size(); i++) {
        int a = edges[i][1], b = edges[i][2];
        if (not ds.merge(a, b)) continue;
        // tree[j][1] = a; tree[j][2] = b; tree[j][0] = edges[i][0];
        cost += edges[i][0];
        j++;
    }
    return (j == n - 1 ? cost : -1);
}