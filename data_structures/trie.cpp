// O(n)
template <int N, char id>
struct trie {
    struct node {
        array<int, N> c;
        int cnt = 0, end = 0;
        node() { c.fill(-1); }
    };
    vector<node> d;

    int add_node() {
        d.push_back(node());
        return d.size() - 1;
    }

    trie() { add_node(); }

    void insert(const string& s) {
        int n = s.size(), rt = 0;
        for (int i = 0; i < n; i++) {
            if (d[rt].c[s[i] - id] == -1) 
                d[rt].c[s[i] - id] = add_node();
            rt = d[rt].c[s[i] - id];
            d[rt].cnt++;
            if (i + 1 == n) d[rt].end++;
        }
    }

    void erase(const string& s) {
        int n = s.size(), rt = 0;
        for (int i = 0; i < n; i++) {
            rt = d[rt].c[s[i] - id];
            d[rt].cnt--;
            if (i + 1 == n) d[rt].end--;
        }
    }
};