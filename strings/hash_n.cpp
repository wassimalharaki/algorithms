mt19937_64 gen(random_device{}());
uniform_int_distribution<int> dist(27, 1e8);

const array<const int, 12> M{
    725676706, 861999952, 544342990,
    116146196, 731624326, 677831782,
    113463194, 459963832, 684889526,
    275296802, 419643634, 299424160
};
vector<int> B;
array<vector<int>, 12> p;

// O(N * n), O(N)
template <int N>
struct hash_n {
    array<vector<int>, N> h;

    hash_n(const string& s) {
        while (B.size() < 12)
            B.push_back(dist(gen));

        for (int i = 0; i < N; i++) {
            if (p[i].empty()) p[i] = {1};
            while (p[i].size() < s.size() + 1)
                p[i].push_back(p[i].back() * B[i] % M[i]);

            h[i].resize(s.size() + 1, 0);
            for (int j = 0; j < (int) s.size(); j++)
                h[i][j + 1] = (h[i][j] * B[i] % M[i] + s[j]) % M[i];
        }
    }

    array<int, N> get_hash_arr(int l, int r) {
        array<int, N> x;
        for (int i = 0; i < N; i++) {
            x[i] = h[i][r] - h[i][l] * p[i][r - l];
            x[i] = (x[i] % M[i] + M[i]) % M[i];
        }
        return x;
    }

    int get_hash_val(int l, int r) {
        int x = 0;
        for (int i = 0; i < N; i++) {
            int y = h[i][r] - h[i][l] * p[i][r - l];
            x ^= (y % M[i] + M[i]) % M[i];
        }
        return x;
    }
};