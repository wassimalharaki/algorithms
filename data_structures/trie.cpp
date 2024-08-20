// O(n)
template <int N, char id>
struct trie {
    struct node {
        array<int, N> c;
        int cnt = 0, end = 0;
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
            if (not d[rt].c[s[i] - id])
                d[rt].c[s[i] - id] = add_node();
            rt = d[rt].c[s[i] - id];
            d[rt].cnt++;
        }
        d[rt].end++;
    }

    void erase(const string& s) {
        int n = s.size(), rt = 0;
        for (int i = 0; i < n; i++) {
            rt = d[rt].c[s[i] - id];
            d[rt].cnt--;
        }
        d[rt].end--;
    }
};