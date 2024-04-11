#include <bits/stdc++.h>
using namespace std;
#define int long long
#define nl '\n'

random_device rd;
mt19937_64 gen(rd());
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
        int n = 100;
        cout << n << nl;
        for (int i = 1; i <= n; i++) {
            cout << random(1, 100) << " ";
        }
        cout << nl;
    }
}