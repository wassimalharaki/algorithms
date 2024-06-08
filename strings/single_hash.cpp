#include <bits/stdc++.h>
using namespace std;
#define int long long

mt19937_64 gen(random_device{}());
const int M = 1e9 + 9;
const int B = uniform_int_distribution<int>(0, M - 1)(gen);
vector<int> p{1};

// O(n)
struct single_hash {
    vector<int> h;

    single_hash(string& s) {
        h.resize(s.size() + 1);
        
        while (p.size() < s.size())
            p.push_back(p.back() * B % M);

        for (int i = 0; i < s.size(); i++)
            h[i + 1] = (h[i] * B % M + s[i]) % M;
    }

    int get_hash(int l, int r) {
        int x = h[r] - h[l] * p[r - l];
        return (x % M + M) % M;
    }
};