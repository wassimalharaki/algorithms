/**
 * O(N^2 * log N)
 */
#include <bits/stdc++.h>
using namespace std;

const double PI = 2 * acos(0.0);

long long sqr(int x) {
    return ((long long) x) * x;
}

// (A[i], B[i]) is the i-th arc erased from the circle
// all arcs should lie in [-PI,PI]
// and should be treated counter-clockwise
vector<double> A, B;

// adding arc to the list of erased arcs
void add(double a1, double a2) {
    // for some sides of rectangle it is possible that a1 > a2
    // in this case we add 2 * PI to a2
    if (a1 > a2) {
        a2 += 2 * PI;
    }

    // arc should lie in (-PI, PI]
    if (a2 > PI) {

        // hence here we divide it in two arcs [a1,PI] and [-PI,a2-2*PI]
        // as arc should be treated counter-clockwise
        A.push_back(a1), B.push_back(PI);
        A.push_back(-PI), B.push_back(a2 - 2 * PI);
    } else if (a1 <= -PI) {

        // and here we divide it in two arcs [a1+2*PI,PI] and [-PI,a2]
        A.push_back(a1 + 2 * PI), B.push_back(PI);
        A.push_back(-PI), B.push_back(a2);
    } else {

        // and here we can add as it is
        A.push_back(a1), B.push_back(a2);
    }
}

int main() {

    // input number of tests
    int T;
    cin >> T;

    // loop over tests
    for (int t = 0; t < T; ++t) {

        // input W, H, N
        int W, H, N;
        cin >> W >> H >> N;

        // scaling W and H
        W *= 100;
        H *= 100;

        // create arrays for X, Y, R
        vector<int> X(N), Y(N), R(N);

        // loop over circles in input
        for(int i = 0; i < N; ++i) {

            // input Xi, Yi, Ri as temporary double variables
            double Xi, Yi ,Ri;
            cin >> Xi >> Yi >> Ri;

            // scaling and rounding them to be integers
            X[i] = int(100 * Xi + 0.5);
            Y[i] = int(100 * Yi + 0.5);
            R[i] = int(100 * Ri + 0.5);
        }

        // the answer to the problem, init by zero
        double res = 0;

        // loop over circles
        for (int i = 0; i < N; ++i) {

            // shortcut for parameters of the i-th circle
            int Xi = X[i], Yi = Y[i], Ri = R[i];

            // circle lies outside the rectangle and do not intersect its sides
            // without this check we will count it to the answer :)
            if (Xi - W >= Ri || Yi - H >= Ri) {
                continue;
            }

            // init vectors A and B
            A.clear(); B.clear();

            // left side
            if (Xi < Ri) {
                // just draw the picture to see why exactly this arc is erased in this case
                // the same should be applied for other sides
                double ang = acos(-double(Xi) / Ri);
                add(ang, -ang);
            }

            // right side
            if (abs(W - Xi) < Ri) {
                // note that we don't need to consider two cases
                // whether center is inside or outside the rectangle
                double ang = acos(double(W - Xi) / Ri);
                add(-ang, ang);
            }

            // down side
            if (Yi < Ri) {
                double ang = asin(double(Yi) / Ri);
                add(ang - PI, -ang);
            }

            // up side
            if (abs(H - Yi) < Ri) {
                double ang = asin(double(H - Yi) / Ri);
                add(ang, PI - ang);
            }

            // needed to check whether the i-th circle is covered completely by some other circle
            bool covered = false;

            // loop over other circle
            for (int j = 0; j < N; ++j) {

                // we don't need to intersect the circle with itself
                if (i == j) {
                    continue;
                }

                // shortcut for parameters of the j-th circle
                int Xj = X[j], Yj = Y[j], Rj = R[j];

                // the square of distance between centers
                long long Dij = sqr(Xj - Xi) + sqr(Yj - Yi);

                // j-th circle covers i-th circle
                if (Rj >= Ri && sqr(Rj - Ri) >= Dij) {
                    covered = true;
                    break;
                }

                // i-th and j-th circles do not interect at two points
                if(Dij >= sqr(Ri + Rj) || Ri >= Rj && sqr(Ri - Rj) >= Dij) {
                    continue;
                }

                // the polar angle of the center of the j-th center
                // with respect to the center of the i-th circle
                double phi = atan2(double(Yj - Yi), double(Xj - Xi));

                // If Oi is the center of the i-th circle
                // Oj is the center of the j-th circle
                // and Aij is one of the intersection points
                // then alp is the angle(Oj,Oi,Aij)
                // we apply cosine theorem to the triangle(Oj,Oi,Aij) to find it
                double alp = acos((sqr(Ri) + Dij - sqr(Rj)) / (2 * Ri * sqrt(double(Dij))));

                // the j-th circle cut exactly this arc from the i-th circle
                add(phi - alp, phi + alp);
            }

            // if circle was completely covered we skip it
            // so actually all found arc are useless now
            if (covered) {
                continue;
            }

            // visible is the length of the visible part of the i-th circle
            // measured in angles, initially we see all circle
            double visible = 2 * PI;

            // then if at least one arc was erased we process arcs
            if (A.size() > 0) {

                // we sort separately left and right end of the arcs
                // it is a trick that based on the fact that the union
                // of all arcs remain the same after sorting
                sort(A.begin(), A.end());
                sort(B.begin(), B.end());

                // (curA, curB) is the current separate segment in the union of arcs
                // it is initialized by the first arc
                double curA = A[0];
                double curB = B[0];

                // now we loop over all remaining arcs
                for (int j = 1; j < A.size(); ++j) {

                    // if right end of the current separate segment
                    // strictly less than left end of the current arc
                    // then this arc belongs to the next separate segment
                    // we don't need to use epsilon stuff here
                    // since it does not change the answer
                    if (curB < A[j]) {

                        // hence we subtract the current segment from the length of visible part
                        visible -= curB - curA;

                        // and set the new left end of the new separate segment
                        curA = A[j];
                    }
                    // the right end of the separate segment is always equals to B[j]
                    // this convenience was brought to us by separate sorting of ends
                    curB = B[j];
                }

                // finally we should subtract the segment in the end since
                // when we finish the cycle this segment remains un-processed
                visible -= curB - curA;
            }

            // we multiply visible by the radius of i-th circle and add to the result
            res += Ri * visible;
        }

        // output the result
        printf("%.7lf\n", res * 0.01);
    }
    return 0;
}
