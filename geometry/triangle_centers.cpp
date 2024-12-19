#include "./point.cpp"

P bary(const P& A, const P& B, const P& C, double a, double b, double c) {
    return (A*a + B*b + C*c) / (a + b + c);
}

P centroid(const P& A, const P& B, const P& C) {
    // geometric center of mass
    return bary(A, B, C, 1, 1, 1);
}

P circumcenter(const P& A, const P& B, const P& C) {
    // intersection of perpendicular bisectors
    double a = B.dist2(C), b = C.dist2(A), c = A.dist2(B);
    return bary(A, B, C, a*(b+c-a), b*(c+a-b), c*(a+b-c));
}

P incenter(const P& A, const P& B, const P& C) {
    // intersection of internal angle bisectors
    return bary(A, B, C, B.dist(C), A.dist(C), A.dist(B));
}

P orthocenter(const P& A, const P& B, const P& C) {
    // intersection of altitudes
    double a = B.dist2(C), b = C.dist2(A), c = A.dist2(B);
    return bary(A, B, C, (a+b-c)*(c+a-b), (b+c-a)*(a+b-c), (c+a-b)*(b+c-a));
}

P excenter(const P& A, const P& B, const P& C) {
    // intersection of two external angle bisectors
    double a = B.dist(C), b = A.dist(C), c = A.dist(B);
    return bary(A, B, C, -a, b, c);

    //// NOTE: there are three excenters
    // return bary(A, B, C, a, -b, c);
    // return bary(A, B, C, a, b, -c);
}
