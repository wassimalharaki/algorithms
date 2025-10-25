template <class E>
struct csr {
    vector<int> start;
    vector<E> elist;
    explicit csr(int n, const vector<pair<int, E>>& edges) {
        elist.resize(edges.size());
        start.resize(n + 1);

        for (auto e : edges)
            start[e.first + 1]++;
        for (int i = 1; i <= n; i++)
            start[i] += start[i - 1];

        auto counter = start;
        for (auto e : edges)
            elist[start[e.first]++] = e.second;
    }
};

// O(F(V + E)log(V + E))
template <class Cap, class Cost>
struct mcf_graph {
    struct edge {
        int from, to;
        Cap cap, flow;
        Cost cost;
    };

    struct _edge {
        int to, rev;
        Cap cap;
        Cost cost;
    };

    int _n;
    vector<edge> _edges;

    mcf_graph() {}
    mcf_graph(int n) : _n(n) {}

    vector<pair<Cap, Cost>> slope(csr<_edge>& g, int s, int t, Cap flow_limit) {
        vector<pair<Cost, Cost>> dual_dist(_n);
        vector<int> prev_e(_n);
        vector<char> vis(_n);
        struct Q {
            Cost key;
            int to;
            bool operator<(Q r) const { return key > r.key; }
        };
        vector<int> que_min;
        vector<Q> que;

        auto dual_ref = [&]() {
            for (int i = 0; i < _n; i++)
                dual_dist[i].second = numeric_limits<Cost>::max();
            fill(vis.begin(), vis.end(), false);
            que_min.clear();
            que.clear();

            size_t heap_r = 0;

            dual_dist[s].second = 0;
            que_min.push_back(s);
            while (not que_min.empty() or not que.empty()) {
                int u;
                if (not que_min.empty()) {
                    u = que_min.back();
                    que_min.pop_back();
                }
                else {
                    while (heap_r < que.size()) {
                        heap_r++;
                        push_heap(que.begin(), que.begin() + heap_r);
                    }
                    u = que.front().to;
                    pop_heap(que.begin(), que.end());
                    que.pop_back();
                    heap_r--;
                }

                if (vis[u]) continue;
                vis[u] = true;
                if (u == t) break;

                Cost dual_u = dual_dist[u].first, dist_u = dual_dist[u].second;
                for (int i = g.start[u]; i < g.start[u + 1]; i++) {
                    auto e = g.elist[i];
                    if (not e.cap) continue;
                    Cost cost = e.cost - dual_dist[e.to].first + dual_u;
                    if (dual_dist[e.to].second - dist_u > cost) {
                        Cost dist_to = dist_u + cost;
                        dual_dist[e.to].second = dist_to;
                        prev_e[e.to] = e.rev;
                        if (dist_to == dist_u)
                            que_min.push_back(e.to);
                        else
                            que.push_back(Q{dist_to, e.to});
                    }
                }
            }
            if (not vis[t]) return false;

            for (int u = 0; u < _n; u++)
                if (vis[u])
                    dual_dist[u].first -=
                        dual_dist[t].second - dual_dist[u].second;
            return true;
        };

        Cap flow = 0;
        Cost cost = 0, prev_cost_per_flow = -1;
        vector<pair<Cap, Cost>> result = {{Cap(0), Cost(0)}};
        while (flow < flow_limit) {
            if (not dual_ref()) break;
            Cap c = flow_limit - flow;

            for (int u = t; u != s; u = g.elist[prev_e[u]].to)
                c = min(c, g.elist[g.elist[prev_e[u]].rev].cap);

            for (int u = t; u != s; u = g.elist[prev_e[u]].to) {
                auto& e = g.elist[prev_e[u]];
                e.cap += c;
                g.elist[e.rev].cap -= c;
            }

            Cost d = -dual_dist[s].first;
            flow += c, cost += c * d;
            if (prev_cost_per_flow == d)
                result.pop_back();
            result.push_back({flow, cost});
            prev_cost_per_flow = d;
        }
        return result;
    }

    int add_edge(int from, int to, Cap cap, Cost cost) {
        _edges.push_back({from, to, cap, 0, cost});
        return _edges.size() - 1;
    }

    edge get_edge(int i) { return _edges[i]; }
    vector<edge> edges() { return _edges; }

    pair<Cap, Cost> flow(int s, int t) {
        return flow(s, t, numeric_limits<Cap>::max());
    }
    pair<Cap, Cost> flow(int s, int t, Cap flow_limit) {
        return slope(s, t, flow_limit).back();
    }
    vector<pair<Cap, Cost>> slope(int s, int t) {
        return slope(s, t, numeric_limits<Cap>::max());
    }
    vector<pair<Cap, Cost>> slope(int s, int t, Cap flow_limit) {
        int m = (int) _edges.size();
        vector<int> edge_idx(m);

        auto g = [&]() {
            vector<int> degree(_n), redge_idx(m);
            vector<pair<int, _edge>> elist;
            elist.reserve(2 * m);
            for (int i = 0; i < m; i++) {
                auto e = _edges[i];
                edge_idx[i] = degree[e.from]++;
                redge_idx[i] = degree[e.to]++;
                elist.push_back({e.from, {e.to, -1, e.cap - e.flow, e.cost}});
                elist.push_back({e.to, {e.from, -1, e.flow, -e.cost}});
            }
            auto _g = csr<_edge>(_n, elist);
            for (int i = 0; i < m; i++) {
                auto e = _edges[i];
                edge_idx[i] += _g.start[e.from];
                redge_idx[i] += _g.start[e.to];
                _g.elist[edge_idx[i]].rev = redge_idx[i];
                _g.elist[redge_idx[i]].rev = edge_idx[i];
            }
            return _g;
        }();

        auto result = slope(g, s, t, flow_limit);
        for (int i = 0; i < m; i++)
            _edges[i].flow = _edges[i].cap - g.elist[edge_idx[i]].cap;
        return result;
    }
};