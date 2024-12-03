/**
 * 16:14:04 11/10/24
 * circle
 */
// ./ICPC/Geometry/OrganizedTemplates/circle.cpp
#include "./polygon.cpp"

// Ps returned in cw order
vector<P> circle_circle_intersection(const P& a, ftype r, const P& b, ftype R) {
    if (a == b and sign(r - R) == 0) return {P(inf, inf)};
    vector<P> ret;
    ftype d = a.dist(b);
    if (ls(r + R, d) or ls(d + min(r,  R), max(r, R))) return ret;
    ftype x = (d * d - R * R + r * r) / (2 * d);
    ftype y = sqrt(r * r - x * x);
    P v = (b - a) / d;
    ret.push_back(a + v * x  + v.rotateccw(pi / 2) * y);
    if (ls(0, y)) ret.push_back(a + v * x - v.rotateccw(pi / 2) * y);
    return ret;
}

vector<P> circle_circle_intersection(const Circle& c0, const Circle& c1) {
    return circle_circle_intersection(c0.F, c0.S, c1.F, c1.S);
}

vector<P> circle_line_intersection(const P& c, ftype r, P a, P b) {
    vector<P> ret;
    b = b - a; a = a - c;
    ftype A = b.dot(b), B = a.dot(b);
    ftype C = a.dot(a) - r * r, D = B * B - A * C;
    if (ls(D, 0)) return ret;
    ret.push_back(c + a + b * (-B + sqrt(D + eps)) / A);
    if (ls(0, D)) ret.push_back(c + a + b * (-B - sqrt(D)) / A);
    return ret;
}

vector<P> circle_seg_intersection(const P& c, ftype r, const P& a, const P& b) {
    auto li = circle_line_intersection(c, r, a, b);
    if (li.empty()) return li;
    V v{a, b};
    vector<P> res;
    while (not li.empty()) {
        if (v.on(li.back())) res.push_back(li.back());
        li.pop_back();
    }
    return res;
}

vector<P> circle_poly_intersection(const Polygon& p, const P& c, ftype r) {
    vec<P> res;
    for (const auto& v: p.segments()) {
        vec<P> inter = circle_seg_intersection(c, r, v.s, v.d);
        for (const P& pnt: inter) if (res.empty() or res.back() != pnt) res.push_back(pnt);
    }
    return res;
}

ftype circle_circle_area(const P& a, ftype r1, const P& b, ftype r2) {
    ftype d = (a - b).norm();
    if (ls(r1 + r2, d)) return 0;
    if (ls(r1 + d, r2)) return pi * r1 * r1;
    if (ls(r2 + d, r1)) return pi * r2 * r2;
    ftype theta_1 = acos((r1 * r1 + d * d - r2 * r2) / (2 * r1 * d)),
          theta_2 = acos((r2 * r2 + d * d - r1 * r1)/(2 * r2 * d));
    return r1 * r1 * (theta_1 - sin(2 * theta_1)/2.) + r2 * r2 * (theta_2 - sin(2 * theta_2)/2.);
}

// tangent lines from point q to the circle
int tangent_lines_from_point(const P& p, ftype r, const P& q, V& u, V& v) {
    int x = sign(p.dist2(q) - r * r);
    if (x < 0) return 0; // point in cricle
    if (x == 0) { // point on circle
        u = V(q, q + (q - p).normal_l());
        v = u;
        return 1;
    }
    ftype d = p.dist(q), l = r * r / d, h = sqrt(r * r - l * l);
    u = V(q, p + ((q - p).truncate(l) + ((q - p).normal_l().truncate(h))));
    v = V(q, p + ((q - p).truncate(l) + ((q - p).normal_r().truncate(h))));
    return 2;
}

// returns outer tangents line of two circles
// if inner == 1 it returns inner tangent lines
int tangents_lines_from_circle(const P& c1, ftype r1, const P& c2, ftype r2, bool inner, V& u, V& v) {
    if (inner) r2 = -r2;
    P d = c2 - c1;
    ftype dr = r1 - r2, d2 = d.norm(), h2 = d2 - dr * dr;
    if (eq(d2, 0) or ls(h2, 0)) {
        assert(not eq(h2, 0));
        return 0;
    }
    vector<V> out;
    for (int tmp: {-1, 1}) {
        P l = (d * dr + d.normal_l() * sqrt(h2) * tmp) / d2;
        out.push_back({c1 + l * r1, c2 + l * r2});
    }
    u = out[0];
    if (out.size() == 2) v = out[1];
    return 1 + ls(0, h2);
}

// find a circle of radius r that contains as many points as possible
// O(n^2 log n);
double maximum_circle_cover(vector<P> p, double r, Circle& c) {
    int n = p.size();
    int ans = 0;
    int id = 0; double th = 0;
    for (int i = 0; i < n; ++i) {
        // maximum circle cover when the circle goes through this point
        vector<pair<double, int>> events = {{-pi, +1}, {pi, -1}};
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            double d = p[i].dist(p[j]);
            if (d > r * 2) continue;
            double dir = (p[j] - p[i]).a();
            double ang = acos(d / 2 / r);
            double st = dir - ang, ed = dir + ang;
            if (st > pi) st -= pi * 2;
            if (st <= -pi) st += pi * 2;
            if (ed > pi) ed -= pi * 2;
            if (ed <= -pi) ed += pi * 2;
            events.push_back({st - eps, +1}); // take care of precisions!
            events.push_back({ed, -1});
            if (st > ed) {
                events.push_back({-pi, +1});
                events.push_back({+pi, -1});
            }
        }
        sort(events.begin(), events.end());
        int cnt = 0;
        for (auto &&e: events) {
            cnt += e.second;
            if (cnt > ans) {
                ans = cnt;
                id = i; th = e.first;
            }
        }
    }
    P w = P(p[id].x + r * cos(th), p[id].y + r * sin(th));
    c = Circle(w, r); //best_circle
    return ans;
}

// a is distance along major axis, b is distance along minor axis, center is (h, k)
// ellipse aligned with coordinate axis.
// To find the solution for rotated ellipses, rotate the whole system and then re-rotate the resulting points in the opposite direction.
vector<P> ellipse_line_intersection(ftype h, ftype k, ftype a, ftype b, const V& v) {
    vector<P> intersections;
    ftype x1 = v.s.x, y1 = v.s.y, x2 = v.d.x, y2 = v.d.y;

    // Calculate line equation parameters (slope m and intercept c)
    ftype dx = x2 - x1;
    ftype dy = y2 - y1;
    if (eq(dx, 0)) {  // Handle vertical line case
        ftype X = x1 - h;
        ftype A = 1.0 / (b * b);
        ftype B = -2.0 * k / (b * b);
        ftype C = (X * X) / (a * a) + (k * k) / (b * b) - 1.0;

        auto ySolutions = quadratic(A, B, C);
        if (ySolutions) {
            for (ftype Y : ySolutions.value()) {
                intersections.emplace_back(x1, Y + k);
            }
        }
        return intersections;
    }

    // Line is not vertical; calculate slope (m) and intercept (c)
    ftype m = dy / dx;
    ftype c = y1 - m * x1;

    ftype A = 1.0 / (a * a) + (m * m) / (b * b);
    ftype B = (2 * m * (c - k) - 2 * h) / (a * a);
    ftype C = ((c - k) * (c - k)) / (b * b) + (h * h) / (a * a) - 1.0;

    auto xSolutions = quadratic(A, B, C);
    if (xSolutions) {
        for (ftype X : xSolutions.value()) {
            ftype Y = m * X + c;
            intersections.emplace_back(X, Y);
        }
    }

    return intersections;
}

//5 - outside and do not intersect
//4 - intersect outside in one point
//3 - intersect in 2 points
//2 - intersect inside in one point
//1 - inside and do not intersect
int circle_circle_relation(const P& a, double r, const P& b, double R) {
    double d = a.dist(b);
    if (sign(d - r - R) > 0)  return 5;
    if (sign(d - r - R) == 0) return 4;
    double l = fabs(r - R);
    if (sign(d - r - R) < 0 && sign(d - l) > 0) return 3;
    if (sign(d - l) == 0) return 2;
    if (sign(d - l) < 0) return 1;
    assert(0); return -1;
}

// The case for ellipse-ellipse intersection uses newton's method
// The following code is JS
//
//var MAX_ITERATIONS = 10
//var innerPolygonCoef, outerPolygonCoef, initialized
//
//function initialize()
//{
//    innerPolygonCoef = []
//    outerPolygonCoef = []
//    for (var t = 0; t <= MAX_ITERATIONS; t++)
//    {
//        var numNodes = 4 << t
//        innerPolygonCoef[t] = 0.5 / Math.cos(4 * Math.acos(0) / numNodes)
//        outerPolygonCoef[t] = 0.5 / (Math.cos(2 * Math.acos(0) / numNodes) * Math.cos(2 * Math.acos(0) / numNodes))
//    }
//    initialized = true
//}
//
//function iterate(x, y, c0x, c0y, c2x, c2y, rr)
//{
//    for (var t = 1; t <= MAX_ITERATIONS; t++)
//    {
//        var c1x = (c0x + c2x) * innerPolygonCoef[t]
//        var c1y = (c0y + c2y) * innerPolygonCoef[t]
//        var tx = x - c1x
//        var ty = y - c1y
//        if (tx * tx + ty * ty <= rr)
//        {
//            return true
//        }
//        var t2x = c2x - c1x
//        var t2y = c2y - c1y
//        if (tx * t2x + ty * t2y >= 0 && tx * t2x + ty * t2y <= t2x * t2x + t2y * t2y &&
//            (ty * t2x - tx * t2y >= 0 || rr * (t2x * t2x + t2y * t2y) >= (ty * t2x - tx * t2y) * (ty * t2x - tx * t2y)))
//        {
//            return true
//        }
//        var t0x = c0x - c1x
//        var t0y = c0y - c1y
//        if (tx * t0x + ty * t0y >= 0 && tx * t0x + ty * t0y <= t0x * t0x + t0y * t0y &&
//            (ty * t0x - tx * t0y <= 0 || rr * (t0x * t0x + t0y * t0y) >= (ty * t0x - tx * t0y) * (ty * t0x - tx * t0y)))
//        {
//            return true
//        }
//        var c3x = (c0x + c1x) * outerPolygonCoef[t]
//        var c3y = (c0y + c1y) * outerPolygonCoef[t]
//        if ((c3x - x) * (c3x - x) + (c3y - y) * (c3y - y) < rr)
//        {
//            c2x = c1x
//            c2y = c1y
//            continue
//        }
//        var c4x = c1x - c3x + c1x
//        var c4y = c1y - c3y + c1y
//        if ((c4x - x) * (c4x - x) + (c4y - y) * (c4y - y) < rr)
//        {
//            c0x = c1x
//            c0y = c1y
//            continue
//        }
//        var t3x = c3x - c1x
//        var t3y = c3y - c1y
//        if (ty * t3x - tx * t3y <= 0 || rr * (t3x * t3x + t3y * t3y) > (ty * t3x - tx * t3y) * (ty * t3x - tx * t3y))
//        {
//            if (tx * t3x + ty * t3y > 0)
//            {
//                if (Math.abs(tx * t3x + ty * t3y) <= t3x * t3x + t3y * t3y || (x - c3x) * (c0x - c3x) + (y - c3y) * (c0y - c3y) >= 0)
//                {
//                    c2x = c1x
//                    c2y = c1y
//                    continue
//                }
//            } else if (-(tx * t3x + ty * t3y) <= t3x * t3x + t3y * t3y || (x - c4x) * (c2x - c4x) + (y - c4y) * (c2y - c4y) >= 0)
//            {
//                c0x = c1x
//                c0y = c1y
//                continue
//            }
//        }
//        return false
//    }
//    return false // Out of iterations so it is unsure if there was a collision. But have to return something.
//}
//
//// Test for collision between an ellipse of horizontal radius w0 and vertical radius h0 at (x0, y0) and
//// an ellipse of horizontal radius w1 and vertical radius h1 at (x1, y1)
//function ellipseEllipse(x0, y0, w0, h0, x1, y1, w1, h1)
//{
//    if (!initialized)
//    {
//        initialize()
//    }
//
//    var x = Math.abs(x1 - x0) * h1
//    var y = Math.abs(y1 - y0) * w1
//    w0 *= h1
//    h0 *= w1
//    var r = w1 * h1
//
//    if (x * x + (h0 - y) * (h0 - y) <= r * r || (w0 - x) * (w0 - x) + y * y <= r * r || x * h0 + y * w0 <= w0 * h0
//        || ((x * h0 + y * w0 - w0 * h0) * (x * h0 + y * w0 - w0 * h0) <= r * r * (w0 * w0 + h0 * h0) && x * w0 - y * h0 >= -h0 * h0 && x * w0 - y * h0 <= w0 * w0))
//    {
//        return true
//    }
//    else
//    {
//        if ((x - w0) * (x - w0) + (y - h0) * (y - h0) <= r * r || (x <= w0 && y - r <= h0) || (y <= h0 && x - r <= w0))
//        {
//            return iterate(x, y, w0, 0, 0, h0, r * r)
//        }
//        return false
//    }
//}
//
//// Test for collision between an ellipse of horizontal radius w and vertical radius h at (x0, y0) and
//// a circle of radius r at (x1, y1)
//function ellipseCircle(x0, y0, w, h, x1, y1, r)
//{
//    if (!initialized)
//    {
//        initialize()
//    }
//    var x = Math.abs(x1 - x0)
//    var y = Math.abs(y1 - y0)
//
//    if (x * x + (h - y) * (h - y) <= r * r || (w - x) * (w - x) + y * y <= r * r || x * h + y * w <= w * h
//        || ((x * h + y * w - w * h) * (x * h + y * w - w * h) <= r * r * (w * w + h * h) && x * w - y * h >= -h * h && x * w - y * h <= w * w))
//    {
//        return true
//    }
//    else
//    {
//        if ((x - w) * (x - w) + (y - h) * (y - h) <= r * r || (x <= w && y - r <= h) || (y <= h && x - r <= w))
//        {
//            return iterate(x, y, w, 0, 0, h, r * r)
//        }
//        return false
//    }
//}
