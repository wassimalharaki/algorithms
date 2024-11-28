/**
 * 12:41:52 11/24/24
 * slope_trick
 */
// ./algorithms/geometry/OrganizedTemplates/slope_trick.cpp
#include <bits/stdc++.h>
using namespace std;

// https://codeforces.com/blog/entry/47821
void solve() {
    int n, t;
    long long ans = 0;
    std::priority_queue<int> Q;
    scanf("%d%d", &n, &t);
    Q.push(t);
    for(int i=1; i<n; i++)
    {
        scanf("%d", &t); t-=i;
        Q.push(t);
        if(Q.top() > t)
        {
            ans += Q.top() - t;
            Q.pop();
            Q.push(t);
        }
    }
    printf("%lld", ans);
}
