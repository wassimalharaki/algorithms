#include <bits/stdc++.h>
using namespace std;
#define int long long
#define nl '\n'
#define v vector

mt19937_64 gen(random_device{}());
uniform_int_distribution<int> dist(1, 1);

int random(int a, int b) {
    dist.param(uniform_int_distribution<int>::param_type(a, b));
    return dist(gen);
}

signed main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);

    freopen("input.txt", "w", stdout);
    int T = 1000;
    cout << T << nl;
    for (int t = 1; t <= T; t++) {
        int n = random(2, 10), m = random(1, n * (n - 1) / 2);
        cout << n << " " << m << nl;
        for (int i = 0; i < m; i++) {
            cout << random(1, n) << " " << random(1, n) << nl;
        }
    }
}