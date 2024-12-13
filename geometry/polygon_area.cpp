struct point {
    int x, y;
};

// must be cyclically sorted
// O(n)
double area(const vector<point>& pts) {
    double res = 0;
    for (int i = 0; i < pts.size(); i++) {
        point p = i ? pts[i - 1] : pts.back();
        point q = pts[i];
        res += (p.x - q.x) * (p.y + q.y);
    }
    return fabs(res) / 2;
}