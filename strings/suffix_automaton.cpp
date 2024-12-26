// O(nN)
template<int N, char id>
struct suffix_automaton {
    struct node {
        int sz = 0, link = -1;
        array<int, N> go;
        node() { go.fill(-1); }
    };
    vector<node> d;
    int last = 0;

    int add_node() {
        d.push_back(node());
        return d.size() - 1;
    }

    suffix_automaton() { add_node(); }

    suffix_automaton(const string& s) {
        add_node();
        for (const char& c : s)
            extend(c);
    }

    void extend(char c) {
        int curr = add_node();
        d[curr].sz = d[last].sz + 1;

        int p = last;
        while (p != -1 and d[p].go[c - id] == -1) {
            d[p].go[c - id] = curr;
            p = d[p].link;
        }

        if (p == -1)
            d[curr].link = 0;
        else {
            int q = d[p].go[c - id];
            if (d[p].sz + 1 == d[q].sz)
                d[curr].link = q;
            else {
                int clone = add_node();
                d[clone].sz = d[p].sz + 1;
                d[clone].go = d[q].go;
                d[clone].link = d[q].link;
                while (p != -1 and d[p].go[c - id] == q) {
                    d[p].go[c - id] = clone;
                    p = d[p].link;
                }
                d[q].link = d[curr].link = clone;
            }
        }
        last = curr;
    }
};