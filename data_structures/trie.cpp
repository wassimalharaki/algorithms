// O(n * N)
template <int N, char id>
struct trie {
    struct node {
        array<int, N> go;
        int cnt = 0, end = 0;
    };
    vector<node> d;

    int add_node() {
        d.push_back(node());
        return d.size() - 1;
    }

    trie() { add_node(); }

    void insert(const string& s) {
        int rt = 0;
        for (const char& c : s) {
            if (not d[rt].go[c - id])
                d[rt].go[c - id] = add_node();
            rt = d[rt].go[c - id];
            d[rt].cnt++;
        }
        d[rt].end++;
    }

    void erase(const string& s) {
        int rt = 0;
        for (const char& c : s) {
            rt = d[rt].go[c - id];
            d[rt].cnt--;
        }
        d[rt].end--;
    }
};