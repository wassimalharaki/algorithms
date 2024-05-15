#include <bits/stdc++.h>
using namespace std;

template <class Cap> struct mf_graph {
    int n;
    struct _edge {
        int to, rev;
        Cap cap;
    };
    vector<pair<int, int>> pos;
    vector<vector<_edge>> g;
    
    mf_graph(int n) : n(n), g(n) {}

    int add_edge(int from, int to, Cap cap) {
        int m = pos.size();
        pos.push_back({from, g[from].size()});
        int from_id = g[from].size();
        int to_id = g[to].size();
        if (from == to) to_id++;
        g[from].push_back(_edge{to, to_id, cap});
        g[to].push_back(_edge{from, from_id, 0});
        return m;
    }

    struct edge {
        int from, to;
        Cap cap, flow;
    };

    edge get_edge(int i) {
        int m = pos.size();
        auto e = g[pos[i].first][pos[i].second];
        auto re = g[e.to][e.rev];
        return edge{pos[i].first, e.to, e.cap + re.cap, re.cap};
    }

    vector<edge> edges() {
        int m = pos.size();
        vector<edge> result;
        for (int i = 0; i < m; i++)
            result.push_back(get_edge(i));
        return result;
    }

    void change_edge(int i, Cap new_cap, Cap new_flow) {
        int m = pos.size();
        auto& e = g[pos[i].first][pos[i].second];
        auto& re = g[e.to][e.rev];
        e.cap = new_cap - new_flow;
        re.cap = new_flow;
    }

    Cap flow(int s, int t) {
        return flow(s, t, numeric_limits<Cap>::max());
    }

    Cap flow(int s, int t, Cap flow_limit) {
        vector<int> level(n), iter(n);
        queue<int> q;

        auto bfs = [&]() {
            fill(level.begin(), level.end(), -1);
            level[s] = 0;
            while (q.size()) q.pop();
            q.push(s);
            while (!q.empty()) {
                int v = q.front();
                q.pop();
                for (auto e : g[v]) {
                    if (e.cap == 0 || level[e.to] >= 0) continue;
                    level[e.to] = level[v] + 1;
                    if (e.to == t) return;
                    q.push(e.to);
                }
            }
        };
        auto dfs = [&](int v, Cap up, auto&& dfs) {
            if (v == s) return up;
            Cap res = 0;
            int level_v = level[v];
            for (int& i = iter[v]; i < g[v].size(); i++) {
                _edge& e = g[v][i];
                if (level_v <= level[e.to] || g[e.to][e.rev].cap == 0) continue;
                Cap d = dfs(e.to, min(up - res, g[e.to][e.rev].cap), dfs);
                if (d <= 0) continue;
                g[v][i].cap += d;
                g[e.to][e.rev].cap -= d;
                res += d;
                if (res == up) return res;
            }
            level[v] = n;
            return res;
        };

        Cap flow = 0;
        while (flow < flow_limit) {
            bfs();
            if (level[t] == -1) break;
            fill(iter.begin(), iter.end(), 0);
            Cap f = dfs(t, flow_limit - flow, dfs);
            if (!f) break;
            flow += f;
        }
        return flow;
    }

    vector<bool> min_cut(int s) {
        vector<bool> visited(n);
        queue<int> q;
        q.push(s);
        while (!q.empty()) {
            int p = q.front();
            q.pop();
            visited[p] = true;
            for (auto e : g[p]) {
                if (e.cap && !visited[e.to]) {
                    visited[e.to] = true;
                    q.push(e.to);
                }
            }
        }
        return visited;
    }
};