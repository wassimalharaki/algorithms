#include <bits/stdc++.h>
using namespace std;
#define v vector

v<int> comp, comps;
void tarjan(v<v<int>>& adj) {
    int n = adj.size(), curr = 0;
    comp.assign(n, -1);
    v<int> disc(n), vis;

    auto dfs = [&](int u, auto&& dfs) -> int {
        int low = disc[u] = ++curr;
		vis.push_back(u);

		for (int& i : adj[u])
			if (comp[i] == -1)
                low = min(low, disc[i] ?: dfs(i, dfs));

		if (low == disc[u]) {
			comps.push_back(u);
			for (int i = -1; i != u;) {
                i = vis.back();
                comp[i] = u;
                vis.pop_back();
            }
		}
		return low;
    };

    for (int i = 0; i < n; i++)
        if (!disc[i]) dfs(i, dfs);
    reverse(comps.begin(), comps.end());
}