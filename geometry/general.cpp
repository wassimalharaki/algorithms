struct star {
    int n;    // number of sides of the star
    double r; // radius of the circumcircle
    star(int _n, double _r) {
        n = _n;
        r = _r;
    }

    double area() {
        double theta = PI / n;
        double s = 2 * r * sin(theta);
        double R = 0.5 * s / tan(theta);
        double a = 0.5 * n * s * R;
        double a2 = 0.25 * s * s / tan(1.5 * theta);
        return a - n * a2;
    }
};

// given a list of lengths of the sides of a polygon in counterclockwise order
// returns the maximum area of a non-degenerate polygon that can be formed using those lengths
double get_maximum_polygon_area_for_given_lengths(vector<double> v) {
    if (v.size() < 3) {
        return 0;
    }
    int m = 0;
    double sum = 0;
    for (int i = 0; i < v.size(); i++) {
        if (v[i] > v[m]) {
            m = i;
        }
        sum += v[i];
    }
    if (sign(v[m] - (sum - v[m])) >= 0) {
        return 0; // no non-degenerate polygon is possible
    }
    // the polygon should be a circular polygon
    // that is all points are on the circumference of a circle
    double l = v[m] / 2, r = 1e6; // fix it correctly
    int it = 60;
    auto ang = [](double x, double r) { // x = length of the chord, r = radius of the circle
        return 2 * asin((x / 2) / r);
    };
    auto calc = [=](double r) {
        double sum = 0;
        for (auto x: v) {
            sum += ang(x, r);
        }
        return sum;
    };
    // compute the radius of the circle
    while (it--) {
        double mid = (l + r) / 2;
        if (calc(mid) <= 2 * PI) {
            r = mid;
        }
        else {
            l = mid;
        }
    }

    if (calc(r) <= 2 * PI - eps) { // the center of the circle is outside the polygon
        auto calc2 = [&](double r) {
            double sum = 0;
            for (int i = 0; i < v.size(); i++) {
                double x = v[i];
                double th = ang(x, r);
                if (i != m) {
                    sum += th;
                }
                else {
                    sum += 2 * PI - th;
                }
            }
            return sum;
        };
        l = v[m] / 2; r = 1e6;
        it = 60;
        while (it--) {
            double mid = (l + r) / 2;
            if (calc2(mid) > 2 * PI) {
                r = mid;
            }
            else {
                l = mid;
            }
        }
        auto get_area = [=](double r) {
            double ans = 0;
            for (int i = 0; i < v.size(); i++) {
                double x = v[i];
                double area = r * r * sin(ang(x, r)) / 2;
                if (i != m) {
                    ans += area;
                }
                else {
                    ans -= area;
                }
            }
            return ans;
        };
        return get_area(r);
    }
    else { // the center of the circle is inside the polygon
        auto get_area = [=](double r) {
            double ans = 0;
            for (auto x: v) {
                ans += r * r * sin(ang(x, r)) / 2;
            }
            return ans;
        };
        return get_area(r);
    }
}
