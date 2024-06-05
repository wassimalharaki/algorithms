#include <bits/stdc++.h>
using namespace std;

template <class F, class S, S (*op)(S, S), S (*e)()>
struct sparse_segtree {
    struct node {
        F l, r;
        S x = e();
        node *left_child = 0, *right_child = 0, *parent;

        node(F _l, F _r, node *_parent) {
            l = _l, r = _r, parent = _parent;
        }

        node* get_left_child() {
            if (left_child == 0)
                left_child = new node(l, l + r >> 1, this);
            return left_child;
        }

        node* get_right_child() {
            if (right_child == 0)
                right_child = new node(l + r >> 1, r, this);
            return right_child;
        }

        void update() {
            x = op(
                left_child ? left_child -> x : e(),
                right_child ? right_child -> x : e()
            );
        }
    };
    typedef node* nodeptr;
    nodeptr root;

    sparse_segtree(F n) {
        root = new node((F) 0, n, 0);
    }

    S get(F i) {
        nodeptr curr = root;
        while (curr and (i != curr -> l or curr -> l + 1 != curr -> r))
            curr = i < curr -> l + curr -> r >> 1 ?
                curr -> left_child :
                curr -> right_child;
        return curr ? curr -> x : e();
    }

    void set(F i, S x) {
        nodeptr curr = root;
        while (curr -> l != i or curr -> l + 1 != curr -> r)
            curr = i < curr -> l + curr -> r >> 1 ?
                curr -> get_left_child() :
                curr -> get_right_child();
        
        curr -> x = x;
        while (curr -> parent) {
            curr = curr -> parent;
            curr -> update();
        }
    }

    void apply(F i, S x) {
        nodeptr curr = root;
        curr -> x = op(curr -> x, x);
        while (curr -> l != i or curr -> l + 1 != curr -> r) {
            curr = i < curr -> l + curr -> r >> 1 ?
                curr -> get_left_child() :
                curr -> get_right_child();
            curr -> x = op(curr -> x, x);
        }
    }

    S prod(F l, F r) {
        return prod_util(l, r, root);
    }

    S prod_util(F l, F r, nodeptr curr) {
        if (curr == 0 or curr -> r <= l or curr -> l >= r)
            return e();
        if (l <= curr -> l and curr -> r <= r)
            return curr -> x;
        return op(
            prod_util(l, r, curr -> left_child),
            prod_util(l, r, curr -> right_child)
        );
    }
};