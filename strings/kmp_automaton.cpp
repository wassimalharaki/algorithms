// O(nN)
template <int N, char id>
vector<vector<int>> kmp_automaton(string s) {
    s.push_back('#');
    int n = s.size();

    vector<int> pi(n);
    for (int i = 1; i < n; i++) {
        int j = pi[i - 1];
        while (j > 0 and s[i] != s[j])
            j = pi[j - 1];
        j += s[i] == s[j];
        pi[i] = j;
    }

    vector adj(n, vector<int>(N));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < N; j++)
            if (i > 0 and id + j != s[i])
                adj[i][j] = adj[pi[i-1]][j];
            else
                adj[i][j] = i + (id + j == s[i]);
    return adj;
}