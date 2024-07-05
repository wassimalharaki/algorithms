// O(V^3)
vector<vector<int>> floyd_warshall(vector<vector<int>> adj) {
    int n = adj.size();
    for (int k = 0; k < n; k++)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (adj[i][j] > adj[i][k] + adj[k][j])
                    adj[i][j] = adj[i][k] + adj[k][j];
    return adj;
}