// O(log(n))
template <int N>
struct binary_trie {
    struct node {
        array<int, 2> c;
        int cnt = 0;
        node() { c.fill(-1); }
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
            if (d[rt].c[b] == -1)
                d[rt].c[b] = add_node();
            rt = d[rt].c[b];
            d[rt].cnt++;
        }
    }

    void erase(int x) {
        int rt = 0;
        for (int i = N; i >= 0; i--) {
            bool b = (1ll << i) & x;
            rt = d[rt].c[b];
            d[rt].cnt--;
        }
    }
 
    int min_xor(int x) {
        int rt = 0, ans = 0;
        for (int i = N; i >= 0; i--) {
            bool b = (1ll << i) & x;
            if (d[rt].c[b] != -1 and d[d[rt].c[b]].cnt)
                rt = d[rt].c[b];
            else
                ans += 1ll << i, rt = d[rt].c[!b];
        }
        return ans;
    }

    int max_xor(int x) {
        int rt = 0, ans = 0;
        for (int i = N; i >= 0; i--) {
            bool b = !((1ll << i) & x);
            if (d[rt].c[b] != -1 and d[d[rt].c[b]].cnt)
                ans += 1ll << i, rt = d[rt].c[b];
            else
                rt = d[rt].c[!b];
        }
        return ans;
    }
};