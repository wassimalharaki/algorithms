// O(N)
template <int N>
struct binary_trie {
    struct node {
        array<int, 2> go;
        int cnt = 0;
    };
    vector<node> d;

    int add_node() {
        d.push_back(node());
        return d.size() - 1;
    }

    binary_trie() { add_node(); }

    void insert(int x) {
        int rt = 0;
        for (int i = N; i >= 0; i--) {
            bool b = (1ll << i) & x;
            if (not d[rt].go[b])
                d[rt].go[b] = add_node();
            rt = d[rt].go[b];
            d[rt].cnt++;
        }
    }

    void erase(int x) {
        int rt = 0;
        for (int i = N; i >= 0; i--) {
            bool b = (1ll << i) & x;
            rt = d[rt].go[b];
            d[rt].cnt--;
        }
    }

    int min_xor(int x) {
        int rt = 0, ans = 0;
        for (int i = N; i >= 0; i--) {
            bool b = (1ll << i) & x;
            if (d[rt].go[b] and d[d[rt].go[b]].cnt)
                rt = d[rt].go[b];
            else
                ans += 1ll << i, rt = d[rt].go[!b];
        }
        return ans;
    }

    int max_xor(int x) {
        int rt = 0, ans = 0;
        for (int i = N; i >= 0; i--) {
            bool b = !((1ll << i) & x);
            if (d[rt].go[b] and d[d[rt].go[b]].cnt)
                ans += 1ll << i, rt = d[rt].go[b];
            else
                rt = d[rt].go[!b];
        }
        return ans;
    }
};