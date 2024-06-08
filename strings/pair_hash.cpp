#include <bits/stdc++.h>
using namespace std;
#define int long long

mt19937_64 gen(random_device{}());
uniform_int_distribution<int> dist(0, 1e9 + 7);
const array<const int, 2> M{(int) 1e9 + 7, (int) 1e9 + 9};
const array<const int, 2> B{dist(gen), dist(gen)};
vector<vector<int>> p{{1}, {1}};

// O(n)
struct pair_hash {
    vector<vector<int>> h{{}, {}};

    pair_hash(string& s) {
        for (int j = 0; j < 2; j++) {
            h[j].resize(s.size() + 1);
            while (p[j].size() < s.size())
                p[j].push_back(p[j].back() * B[j] % M[j]);
        }

        for (int i = 0; i < s.size(); i++)
            for (int j = 0; j < 2; j++)
                h[j][i + 1] = (h[j][i] * B[j] % M[j] + s[i]) % M[j];
    }

    array<int, 2> get_hash(int l, int r) {
        int x = h[0][r] - h[0][l] * p[0][r - l];
        int y = h[1][r] - h[1][l] * p[1][r - l];
        return {(x % M[0] + M[0]) % M[0], (y % M[1] + M[1]) % M[1]};
    }
};