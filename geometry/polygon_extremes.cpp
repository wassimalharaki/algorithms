#include "./polygon.cpp"

using namespace std;
template <class Function>
int extremeVertex(const Polygon& poly, Function direction) {
    int n = static_cast<int>(poly.ps.size()), left = 0, leftSgn;
    auto vertexCmp = [&poly, direction](int i, int j) {
        return sign(ccw(direction(poly.ps[j]), poly.ps[j] - poly.ps[i])); };
    auto isExtreme = [n, vertexCmp](int i, int& iSgn) {
        return (iSgn = vertexCmp((i + 1) % n, i)) >= 0 && vertexCmp(i, i ? i - 1 : (n - 1)) < 0; };
    for (int right = isExtreme(0, leftSgn) ? 1 : n; left + 1 < right;) {
        int middle = (left + right) / 2, middleSgn;
        if (isExtreme(middle, middleSgn)) return middle;
        if (leftSgn != middleSgn ? leftSgn < middleSgn
                                 : leftSgn == vertexCmp(left, middle)) right = middle;
        else left = middle, leftSgn = middleSgn;
    }
    return left;
}

pair<int, int> tangentsConvex(const P& point, const Polygon& poly) {
    return {
            extremeVertex(poly, [&point](const P& q) { return q - point; }),
            extremeVertex(poly, [&point](const P& q) { return point - q; })
    };
}

ftype maxDist2(const Polygon& poly) {
    int n = static_cast<int>(poly.ps.size());
    ftype res = 0;
    for (int i = 0, j = n < 2 ? 0 : 1; i < j; ++i)
        for (;; j = (j + 1) % n) {
            res = max(res, poly.ps[i].dist2(poly.ps[j]));
            if (lse(0, ccw(poly.ps[i+1] - poly.ps[i], poly.ps[(j + 1) % n] - poly.ps[j]))) break;
        }
    return res;
}
