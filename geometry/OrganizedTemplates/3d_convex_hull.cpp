#include "./point_3d.cpp"

// https://codeforces.com/blog/entry/81768

// O(n ^ 2)
// A face is represented by the indices of its three points a, b, c.
// It also stores an outward-facing normal vector q
struct face {
    int a, b, c;
    P q;
};

vector<face> hull3(const vector<P> &p) {
    int n = (int)p.size();
    assert(n >= 3);
    vector<face> f;

    // Consider an edge (a->b) dead if it is not a CCW edge of some current face
    // If an edge is alive but not its reverse, this is an exposed edge.
    // We should add new faces on the exposed edges.
    vector<vector<bool>> dead(n, vector<bool>(n, true));
    auto add_face = [&](int a, int b, int c) {
        f.push_back({a, b, c, (p[b] - p[a]) * (p[c] - p[a])});
        dead[a][b] = dead[b][c] = dead[c][a] = false;
    };

    // Initialize the convex hull of the first 3 points as a
    // triangular disk with two faces of opposite orientation
    add_face(0, 1, 2);
    add_face(0, 2, 1);
    for (int i = 3; i < n; i++) {
        // f2 will be the list of faces invisible to the added point p[i]
        vector<face> f2;
        for(face &F : f) {
            if(ls(0, (p[i] - p[F.a]).dot(F.q))) {
                // this face is visible to the new point, so mark its edges as dead
                dead[F.a][F.b] = dead[F.b][F.c] = dead[F.c][F.a] = true;
            }else {
                f2.push_back(F);
            }
        }
        // Add a new face for each exposed edge.
        // Only check edges of alive faces for being exposed.
        f.clear();
        for(face &F : f2) {
            int arr[3] = {F.a, F.b, F.c};
            for (int j = 0; j < 3; j++) {
                int a = arr[j], b = arr[(j + 1) % 3];
                if(dead[b][a]) {
                    add_face(b, a, i);
                }
            }
        }
        f.insert(f.end(), all(f2));
    }
    return f;
}
