// O(n)
template <int N, char id>
struct aho_corasick {
    struct node {
        array<int, N> go, old;
        int link = -1, vis = 0;
        vector<int> data;
    };
    vector<node> d;
    int n = 0, mxH = 0;

    int add_node() {
        d.push_back(node());
        return d.size() - 1;
    }

    aho_corasick() { add_node(); }

    void insert(const string& s, int i) {
        n++, mxH = max(mxH, (int) s.size());
        int rt = 0;
        for (const char& c : s) {
            if (not d[rt].go[c - id])
                d[rt].go[c - id] = add_node();
            rt = d[rt].go[c - id];
        }
        d[rt].data.push_back(i);
    }

    void build() {
        queue<int> q;
        q.push(0);

        while (q.size()) {
            int u = q.front();
            d[u].old = d[u].go;
            q.pop();

            for (int i = 0; i < N; i++) {
                int j = d[u].go[i];
                if (not j) continue;

                if (d[u].link == -1)
                    d[j].link = 0;
                else
                    d[j].link = d[d[u].link].go[i];

                q.push(j);
            }

            if (u)
                for (int i = 0; i < N; i++)
                    if (not d[u].go[i]) 
                        d[u].go[i] = d[d[u].link].go[i];
        }
    }

    void find(const string& s) {
        int rt = 0;
        for (const char& c : s) {
            rt = d[rt].go[c - id];
            d[rt].vis++;
        }
    }

    vector<int> solve() {
        vector<vector<int>> height(mxH + 1);
        auto dfs = [&](int u, int h, auto&& dfs) -> void {
            height[h].push_back(u);
            for (int i = 0; i < N; i++)
                if (d[u].old[i])
                    dfs(d[u].old[i], h - 1, dfs);
        };
        dfs(0, mxH, dfs);

        vector<int> ans(n);
        for (auto& x : height)
            for (int& i : x) {
                for (int& j : d[i].data)
                    ans[j] += d[i].vis;
                if (d[i].link != -1)
                    d[d[i].link].vis += d[i].vis;
            }
        return ans;
    }
};