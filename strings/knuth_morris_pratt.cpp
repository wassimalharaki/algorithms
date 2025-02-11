// O(n)
vector<int> kmp(const string& s) {
    int n = s.size();
    vector<int> pi(n);
    for (int i = 1; i < n; i++) {
        int j = pi[i - 1];
        while (j > 0 and s[i] != s[j])
            j = pi[j - 1];
        j += s[i] == s[j];
        pi[i] = j;
    }
    return pi;
}