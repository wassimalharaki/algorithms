// O(n)
template <int N, char id>
struct aho_corasick {
    struct node {
        array<int, N> c;
        int link = -1, vis = 0;
        int cnt = 0, end = 0;
    };
    vector<node> d;

    int add_node() {
        d.push_back(node());
        return d.size() - 1;
    }

    aho_corasick() { add_node(); }

    void insert(const string& s, int j) {
        int n = s.size(), rt = 0;
        for (int i = 0; i < n; i++) {
            if (not d[rt].c[s[i] - id])
                d[rt].c[s[i] - id] = add_node();
            rt = d[rt].c[s[i] - id];
            d[rt].cnt++;
        }
        d[rt].end++;
    }

    void build() {
        queue<int> q;
        q.push(0);

        while (q.size()) {
            int u = q.front();
            q.pop();

            for (int i = 0; i < N; i++) {
                int j = d[u].c[i];
                if (not j) continue;

                if (d[u].link == -1)
                    d[j].link = 0;
                else
                    d[j].link = d[d[u].link].c[i];

                q.push(j); 
            }

            if (u)
                for (int i = 0; i < N; i++)
                    if (not d[u].c[i]) 
                        d[u].c[i] = d[d[u].link].c[i];
        }
    }

    void find(const string& s) {
        int n = s.size(), rt = 0;
        for (int i = 0; i < n; i++) {
            while (rt and not d[rt].c[s[i] - id])
                rt = d[rt].link;
            rt = d[rt].c[s[i] - id];
            for (int j = rt; j and not d[j].vis; j = d[j].link)
                d[j].vis = 1;
        }
    }
};