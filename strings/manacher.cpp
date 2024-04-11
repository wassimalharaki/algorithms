#include <bits/stdc++.h>
using namespace std;

vector<int> manacher(string& t) {
    string s{'$'};
    for(char& c: t) s += string{'#', c};
    s += string{'#', '^'};

    int n = s.size() - 2, l = 1, r = 1;
    vector<int> p(n + 2);
    for(int i = 1; i <= n; i++) {
        p[i] = max(0, min(r - i, p[l + (r - i)]));
        while(s[i - p[i]] == s[i + p[i]])
            p[i]++;
        if(i + p[i] > r)
            l = i - p[i], r = i + p[i];
        p[i]--;
    }
    return vector<int>(begin(p) + 2, end(p) - 2);
}

//check if [l; r] is a palindrome
bool is_palindrome(vector<int>& man, int l, int r) {
    return man[l + r] >= r - l + 1;
}