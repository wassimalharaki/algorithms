/**
 * 20:14:53 11/7/24
 * Polygon
 */
// ./ICPC/Geometry/OrganizedTemplates/Polygon.cpp
#include "./vector.cpp"

namespace Geometry {
    using Circle = pair<P, ftype>;
    using Triangle = tuple<P, P, P>;

    double heron(const P& p1, const P& p2, const P& p3) {
        double a = p1.dist(p2);
        double b = p1.dist(p3);
        double c = p3.dist(p2);
        double s = (a+b+c)/2.0;
        double area = sqrt( s*(s-a)*(s-b)*(s-c) );
        return area;
    }

    bool pointInTriangle(const P& a, const P& b, const P& c, const P& p) {
        ftype s1 = abs(a.cross(b, c));
        ftype s2 = abs(p.cross(a, b)) + abs(p.cross(b, c)) + abs(p.cross(c, a));
        return eq(s1, s2);
    }

    void prepare(vec<P>& ps, vec<P>& seq, P& translation) {
        int n = (int)ps.size(), pos = 0;
        for (int i = 1; i < n; i++) {
            if (ps[i] < ps[pos]) pos = i;
        }
        rotate(ps.begin(), ps.begin() + pos, ps.end());

        n--;
        seq.resize(n);
        for (int i = 0; i < n; i++) seq[i] = ps[i + 1] - ps[0];
        translation = ps[0];
    }

    int orientation(const P& a, const P& b, const P& c) {
        double v = a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y);
        if (ls(v, 0)) return -1; // clockwise
        if (ls(0, v)) return +1; // counter-clockwise
        return 0;
    }

    namespace bentley_ottoman {


    }

    bool cw(const P& a, const P& b, const P& c, bool include_collinear) {
        int o = orientation(a, b, c);
        return o < 0 or (include_collinear and o == 0);
    }

    bool ccw(const P& a, const P& b, const P& c, bool include_collinear) {
        int o = orientation(a, b, c);
        return o > 0 or (include_collinear and o == 0);
    }

    bool is_inside(const Circle& c, const P& p) {
        return lse(p.dist(c.F), c.S);
    }

    // helper for circle_from
    P ccenter(const P& b, const P& c) {
        double B = b.dot(c), C = c.norm(), D = b.cross(c);
        return {
            (c.y * B - b.y * C) / (2 * D),
            (b.x * C - c.x * B) / (2 * D)
        };
    }

    Circle circle_from(const P& a, const P& b) {
        P c = 0.5 * (a + b);
        return {c, 0.5 * a.dist(b)};
    }

    Circle circle_from(const P& a, const P& b, const P& c) {
        P i = ccenter(b - a, c - a);
        i += a;
        return {i, i.dist(a)};
    }

    // Function to check whether a circle
    // encloses the given points
    bool is_valid(const Circle& c, const vec<P>& ps) {
        for (const P& p: ps) if (not is_inside(c, p)) return false;
        return true;
    }

    Circle min_circle_brute(const vec<P>& ps) {
        assert(ps.size() <= 3);
        if (ps.empty()) return {{}, 0};
        if (ps.size() == 1) return {ps[0], 0};
        if (ps.size() == 2) return circle_from(ps[0], ps[1]);

        // To check if MEC can be determined
        // by 2 points only
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                Circle c = circle_from(ps[i], ps[j]);
                if (is_valid(c, ps)) return c;
            }
        }
        return circle_from(ps[0], ps[1], ps[2]);
    }

    // Returns the MEC using Welzl's algorithm
    // Takes a set of input points P and a set R
    // points on the circle boundary.
    // n represents the number of points in P
    // that are not yet processed.
    Circle welzl_helper(vec<P>& ps, vec<P> R, int n) {
        // Base case when all points processed or |R| = 3
        if (n == 0 or R.size() == 3) {
            return min_circle_brute(R);
        }

        // Pick a random point randomly
        int idx = rand() % n;
        P p = ps[idx];

        // Put the picked point at the end of P
        // since it's more efficient than
        // deleting from the middle of the vector
        swap(ps[idx], ps[n - 1]);

        // Get the MEC circle d from the
        // set of points P - {p}
        Circle d = welzl_helper(ps, R, n - 1);

        // If d contains p, return d
        if (is_inside(d, p)) {
            return d;
        }

        // Otherwise, must be on the boundary of the MEC
        R.push_back(p);

        // Return the MEC for P - {p} and R U {p}
        return welzl_helper(ps, R, n - 1);
    }

    Circle welzl(vec<P> p) {
        shuffle(all(p), mt19937(random_device()()));
        return welzl_helper(p, {}, (int)p.size());
    }

    namespace triangulate {
        struct QuadEdge {
            P origin;
            QuadEdge* rot = nullptr;
            QuadEdge* onext = nullptr;
            bool used = false;
            QuadEdge* rev() const {
                return rot->rot;
            }
            QuadEdge* lnext() const {
                return rot->rev()->onext->rot;
            }
            QuadEdge* oprev() const {
                return rot->onext->rot;
            }
            P dest() const {
                return rev()->origin;
            }
        };

        QuadEdge* make_edge(const P& from, const P& to) {
            auto* e1 = new QuadEdge;
            auto* e2 = new QuadEdge;
            auto* e3 = new QuadEdge;
            auto* e4 = new QuadEdge;
            e1->origin = from;
            e2->origin = to;
            e3->origin = e4->origin = inf_pt;
            e1->rot = e3;
            e2->rot = e4;
            e3->rot = e2;
            e4->rot = e1;
            e1->onext = e1;
            e2->onext = e2;
            e3->onext = e4;
            e4->onext = e3;
            return e1;
        }

        void splice(QuadEdge* a, QuadEdge* b) {
            swap(a->onext->rot->onext, b->onext->rot->onext);
            swap(a->onext, b->onext);
        }

        void delete_edge(QuadEdge* e) {
            splice(e, e->oprev());
            splice(e->rev(), e->rev()->oprev());
            delete e->rev()->rot;
            delete e->rev();
            delete e->rot;
            delete e;
        }

        QuadEdge* connect(QuadEdge* a, QuadEdge* b) {
            QuadEdge* e = make_edge(a->dest(), b->origin);
            splice(e, a->lnext());
            splice(e->rev(), b);
            return e;
        }

        bool left_of(const P& p, QuadEdge* e) {
            return ls(0, p.cross(e->origin, e->dest()));
        }

        bool right_of(const P& p, QuadEdge* e) {
            return ls(p.cross(e->origin, e->dest()), 0);
        }

        bool in_circle(const P& a, const P& b, const P& c, const P& d) {
            auto ang = [](const P& l, const P& mid, const P& r) {
                ftype x = mid.dot(l, r);
                ftype y = mid.cross(l, r);
                atype res = atan2(x, y);
                return res;
            };
            atype kek = ang(a, b, c) + ang(c, d, a) - ang(b, c, d) - ang(d, a, b);
            return ls(0, kek);
        }

        pair<QuadEdge*, QuadEdge*> build_tr(int l, int r, vector<P>& p) {
            if (r - l + 1 == 2) {
                QuadEdge* res = make_edge(p[l], p[r]);
                return make_pair(res, res->rev());
            }
            if (r - l + 1 == 3) {
                QuadEdge *a = make_edge(p[l], p[l + 1]), *b = make_edge(p[l + 1], p[r]);
                splice(a->rev(), b);
                int sg = sign(p[l].cross(p[l + 1], p[r]));
                if (sg == 0)
                    return make_pair(a, b->rev());
                QuadEdge* c = connect(b, a);
                if (sg == 1)
                    return make_pair(a, b->rev());
                else
                    return make_pair(c->rev(), c);
            }
            int mid = (l + r) / 2;
            QuadEdge *ldo, *ldi, *rdo, *rdi;
            tie(ldo, ldi) = build_tr(l, mid, p);
            tie(rdi, rdo) = build_tr(mid + 1, r, p);
            while (true) {
                if (left_of(rdi->origin, ldi)) {
                    ldi = ldi->lnext();
                    continue;
                }
                if (right_of(ldi->origin, rdi)) {
                    rdi = rdi->rev()->onext;
                    continue;
                }
                break;
            }
            QuadEdge* basel = connect(rdi->rev(), ldi);
            auto valid = [&basel](QuadEdge* e) { return right_of(e->dest(), basel); };
            if (ldi->origin == ldo->origin)
                ldo = basel->rev();
            if (rdi->origin == rdo->origin)
                rdo = basel;
            while (true) {
                QuadEdge* lcand = basel->rev()->onext;
                if (valid(lcand)) {
                    while (in_circle(basel->dest(), basel->origin, lcand->dest(),
                                     lcand->onext->dest())) {
                        QuadEdge* t = lcand->onext;
                        delete_edge(lcand);
                        lcand = t;
                    }
                }
                QuadEdge* rcand = basel->oprev();
                if (valid(rcand)) {
                    while (in_circle(basel->dest(), basel->origin, rcand->dest(),
                                     rcand->oprev()->dest())) {
                        QuadEdge* t = rcand->oprev();
                        delete_edge(rcand);
                        rcand = t;
                    }
                }
                if (!valid(lcand) && !valid(rcand))
                    break;
                if (!valid(lcand) ||
                    (valid(rcand) && in_circle(lcand->dest(), lcand->origin,
                                               rcand->origin, rcand->dest())))
                    basel = connect(rcand, basel->rev());
                else
                    basel = connect(basel->rev(), lcand->rev());
            }
            return make_pair(ldo, rdo);
        }

        vector<Triangle> delaunay(vector<P> p) {
            sort(all(p));
            auto res = build_tr(0, (int)p.size() - 1, p);
            QuadEdge* e = res.first;
            vector<QuadEdge*> edges = {e};
            while (ls(e->onext->dest().cross(e->dest(), e->origin), 0))
                e = e->onext;
            auto add = [&p, &e, &edges]() {
                QuadEdge* curr = e;
                do {
                    curr->used = true;
                    p.push_back(curr->origin);
                    edges.push_back(curr->rev());
                    curr = curr->lnext();
                } while (curr != e);
            };
            add();
            p.clear();
            int kek = 0;
            while (kek < (int)edges.size()) {
                if (!(e = edges[kek++])->used)
                    add();
            }
            vector<Triangle> ans;
            for (int i = 0; i < (int)p.size(); i += 3) {
                ans.emplace_back(p[i], p[i + 1], p[i + 2]);
            }
            return ans;
        }
    }
    namespace vornoi {
        // check if two vectors are collinear. It might make sense to use a
        // different EPS here, especially if points have integer coordinates
        bool collinear(const P& a, const P& b) {
            return eq(a.cross(b), 0);
        }


        // intersection point of lines ab and cd. Precondition is that they aren't collinear
        P lineline(const P& a, const P& b, const P& c, const P& d) {
            return a + (b - a) * ((c - a).cross(d - c) / (b - a).cross(d - c));
        }

        // circumcircle of points a, b, c. Precondition is that abc is a non-degenerate triangle.
        P circumcenter(const P& a, P b, P c) {
            b = (a + b) * 0.5;
            c = (a + c) * 0.5;
            return lineline(b, b + (b - a).normal_l(), c, c + (c - a).normal_l());
        }

        // x coordinate of sweep-line
        ftype sweepx;

        // an arc on the beacah line is given implicitly by the focus p,
        // the focus q of the following arc, and the position of the sweep-line.
        struct arc {
            mutable P p, q;
            mutable int id = 0, i;
            arc(const P& p, const P& q, int i) : p(p), q(q), i(i) {}

            // get y coordinate of intersection with following arc.
            // don't question my magic formulas
            ftype gety(ftype x) const {
                if(q.y == INF) return INF;
                x += eps;
                P med = (p + q) * 0.5;
                P dir = (p - med).normal_l();
                ftype D = (x - p.x) * (x - q.x);
                return med.y + ((med.x - x) * dir.x + sqrtl(D) * dir.abs()) / dir.y;
            }
            bool operator<(const ftype &y) const {
                return gety(sweepx) < y;
            }
            bool operator<(const arc &o) const {
                return gety(sweepx) < o.gety(sweepx);
            }
        };

        // the beach line will be stored as a multiset of arc objects
        using beach = multiset<arc, less<>>;

        // an event is given by
        //     x: the time of the event
        //     id: If >= 0, it's a point event for index id.
        //         If < 0, it's an ID for a vertex event
        //     it: if a vertex event, the iterator for the arc to be deleted
        struct event {
            ftype x;
            int id;
            beach::iterator it;
            event(ftype x, int id, beach::iterator it) : x(x), id(id), it(it) {}
            bool operator<(const event &e) const {
                return x > e.x;
            }
        };

        struct fortune {
            beach line; // self explanatory
            vector<pair<P, int>> v; // (point, original index)
            priority_queue<event> Q; // priority queue of point and vertex events
            vector<pii> edges; // delaunay edges
            vector<bool> valid; // valid[-id] == true if the vertex event with corresponding id is valid
            int n, ti{}; // number of points, next available vertex ID
            explicit fortune(vector<P> p) {
                n = (int)p.size();
                v.resize(n);
                for (int i = 0; i < n; i++) v[i] = {p[i], i};
                sort(all(v)); // sort points by coordinate, remember original indices for the delaunay edges
            }

            // update the remove event for the arc at position it
            void upd(beach::iterator it) {
                if(it->i == -1) return; // doesn't correspond to a real point
                valid[-it->id] = false; // mark existing remove event as invalid
                auto a = prev(it);
                if(collinear(it->q - it->p, a->p - it->p)) return; // doesn't generate a vertex event
                it->id = --ti; // new vertex event ID
                valid.push_back(true); // label this ID true
                P c = circumcenter(it->p, it->q, a->p);
                ftype x = c.x + (c - it->p).abs();
                // event is generated at time x.
                // make sure it passes the sweep-line, and that the arc truly shrinks to 0
                // ls == b - a > eps
                // sweepx - x < eps
                if(x > sweepx - eps && a->gety(x) + eps > it->gety(x)) {
                    Q.emplace(x, it->id, it);
                }
            }

            // add Delaunay edge
            void add_edge(int i, int j) {
                if(i == -1 || j == -1) return;
                edges.emplace_back(v[i].second, v[j].second);
            }

            // handle a point event
            void add(int i) {
                P p = v[i].first;
                // find arc to split
                auto c = line.lower_bound(p.y);
                // insert new arcs. passing the following iterator gives a slight speed-up
                auto b = line.insert(c, arc(p, c->p, i));
                auto a = line.insert(b, arc(c->p, p, c->i));
                add_edge(i, c->i);
                upd(a); upd(b); upd(c);
            }

            // handle a vertex event
            void remove(beach::iterator it) {
                auto a = prev(it);
                auto b = next(it);
                line.erase(it);
                a->q = b->p;
                add_edge(a->i, b->i);
                upd(a); upd(b);
            }

            // X is a value exceeding all coordinates
            void solve(ftype X = 1e9) {
                // insert two points that will always be in the beach line,
                // to avoid handling edge cases of an arc being first or last
                X *= 3;
                line.insert(arc(P(-X, -X), P(-X, X), -1));
                line.insert(arc(P(-X, X), P(INF, INF), -1));
                // create all point events
                for (int i = 0; i < n; i++) {
                    Q.emplace(v[i].first.x, i, line.end());
                }
                ti = 0;
                valid.assign(1, false);
                while(!Q.empty()) {
                    event e = Q.top(); Q.pop();
                    sweepx = e.x;
                    if(e.id >= 0) {
                        add(e.id);
                    }else if(valid[-e.id]) {
                        remove(e.it);
                    }
                }
            }
        };
    }

    struct Polygon {
        vec<P> ps;

        P& operator[](size_t i) {
            return ps[i];
        }

        const P& operator[](size_t i) const {
            return ps[i];
        }

        explicit Polygon(int n): ps(n) {}
        explicit Polygon(const vec<P>& ps): ps{ps} {}

        Polygon operator-() const {
            vec<P> t{ps};
            for (P& p: t) p = -p;
            return Polygon{t};
        }

        // -1 inside, 0 on, 1 outside. O(n)
        int pointpoly(const P& p) const {
            int n = (int)ps.size(), cnt = 0;
            bool boundary = false;
            for (int i = 0; i < n; i++) {
                int j = (i + 1) % n;
                V sd = {ps[i], ps[j]};

                if (sd.on(p)) {
                    boundary = true;
                }

                if (
                    (lse(ps[i].x, p.x) and ls(p.x, ps[j].x) and p.cross(ps[i], ps[j]) < 0)
                    or (lse(ps[j].x, p.x) and ls(p.x, ps[i].x) and p.cross(ps[j], ps[i]) < 0)
                ) {
                    cnt++;
                }
            }

            if (boundary) {
                return 0;
            } else if (cnt % 2) {
                return -1;
            } else {
                return 1;
            }
        }

        // check if point is in for convex polygon
        // polygon ordered ccw
        // O(log(n))
        bool pinConvex(P p) {
            static bool prepared = false;
            static vec<P> seq;
            static P translation;
            if (not prepared) {
                prepare(ps, seq, translation);
                prepared = true;
            }

            int n = (int)ps.size();
            p -= translation;
            if (not eq(seq[0].cross(p), 0) and sign(seq[0].cross(p)) != sign(seq[0].cross(seq[n - 1]))) return false;
            if (not eq(seq[n - 1].cross(p), 0) and sign(seq[n - 1].cross(p)) != sign(seq[n - 1].cross(seq[0]))) return false;
            if (eq(seq[0].cross(p), 0)) return lse(p.abs(), seq[0].abs());

            int l = 0, r = n - 1;
            while (r - l > 1) {
                int mid = (l + r) / 2;
                int pos = mid;
                if (lse(0, seq[pos].cross(p))) l = mid;
                else r = mid;
            }
            int pos = l;
            return pointInTriangle(seq[pos], seq[pos + 1], P(0, 0), p);
        }

        bool is_convex() const {
            int n = (int)ps.size();
            ftype prev = 0, cur = 0;
            for (int i = 0; i < n; i++) {
                vec<P> tmp {
                    ps[i],
                    ps[(i + 1) % n],
                    ps[(i + 2) % n]
                };
                cur = tmp[0].cross(tmp[1], tmp[2]);
                if (not eq(cur, 0)) {
                    if (ls(cur * prev, 0)) return false;
                    else prev = cur;
                }
            }
            return true;
        }

        // -1 inside, 0 intersects, 1 outside
        int segpoly(const V& v) const {
            for (const V& s: segments()) {
                if (not s.shares(v) and s.intersection(v).t != InterT::no) return 0;
            }
            return pointpoly(v.midpoint()) <= 0 ? -1 : 1;
        }

        vec<V> segments() const {
            vec<V> res(ps.size());
            for (size_t i = 0; i < ps.size(); i++) {
                res[i] = {ps[i], ps[(i + 1) % ps.size()]};
            }
            return res;
        }

        ftype area() const {
            ftype res = 0;
            for (int i = 0; i < ps.size(); i++) {
                res += ps[i].cross(ps[(i + 1) % ps.size()]);
            }
            return abs(res) / 2;
        }

        // Pick's
        // Given a certain lattice polygon with non-zero area.
        // We denote its area by S, the number of points with integer coordinates lying strictly inside the polygon by
        // I and the number of points lying on polygon sides by B.
        // Picks formula states: S = I + B / 2 - 1
        int latArea() const { // returns 2 * S
            pii p = latticeinon();
            int I = p.S, B = p.F;
            return 2 * I + B - 2;
        }

        // works given lattic points, number of points in and on.
        pii latticeinon() const {
            vec<int> x(ps.size()), y(ps.size());
            for (int i = 0; i < ps.size(); i++) {
                x[i] = (int)ps[i].x;
                y[i] = (int)ps[i].y;
            }

            int bounds = 0;
            for (int i = 0; i < ps.size(); i++) {
                if (x[i + 1] == x[i]) {
                    bounds += abs(y[i + 1] - y[i]);
                } else if (y[i + 1] == y[i]) {
                    bounds += abs(x[i + 1] - x[i]);
                } else {
                    bounds += gcd(abs(x[i + 1] - x[i]), abs(y[i + 1] - y[i]));
                }
            }

            return {(2 * area() + 2 - bounds) / 2, bounds};
        }

        // O(n * m)
        ftype distance(const Polygon& p) const {
            ftype res = INF;
            for (int i = 0; i < ps.size(); i++) {
                V v{ps[i], ps[(i + 1) % ps.size()]};
                for (const P& pnt: p.ps) res = min(res, v.dist(pnt));
            }
            return res;
        }

        ftype distance(const P& q) const {
            ftype r = INF;
            for (const V& v: segments()) r = min(r, v.dist(q));
            return r;
        }

        vec<Triangle> triangulate() const {
            // O(nlog(n)) convex polygon
            using namespace triangulate;
            return delaunay(ps);
        }

        // Note that if Minkowski subtraction of two convex polygons covers (0, 0), then they intersect
        ftype fastDistance(const Polygon& p) const {
            return minkowskiSum(-p).distance({0, 0});
        }

        bool intersect(const Polygon& q) const {
            return minkowskiSum(-q).pointpoly({0, 0}) <= 0;
        }

        // for minkowski
        void reorder_polygon() {
            int pos = 0;
            for (int i = 1; i < (int)ps.size(); i++) {
                if (ls(ps[i].y, ps[pos].y) or (eq(ps[i].y, ps[pos].y) and ls(ps[i].x, ps[pos].x))) pos = i;
            }
            rotate(ps.begin(), ps.begin() + pos, ps.end());
        }

        Polygon minkowskiSum(Polygon q) const {
            Polygon p = *this;
            p.reorder_polygon();
            q.reorder_polygon();
            p.ps.push_back(p.ps[0]);
            p.ps.push_back(p.ps[1]);
            q.ps.push_back(q.ps[0]);
            q.ps.push_back(q.ps[1]);
            vec<P> res;
            int i = 0, j = 0;
            while (i + 2 < (int)p.ps.size() or j + 2 < q.ps.size()) {
                res.push_back(p.ps[i] + q.ps[j]);
                ftype cross = (p.ps[i + 1] - p.ps[i]).cross(q.ps[j + 1] - q.ps[j]);
                if (lse(0, cross) and i + 2 < (int)p.ps.size()) i++;
                if (lse(0, cross) and j + 2 < (int)q.ps.size()) j++;
            }
            return Polygon{res};
        }

        // Minimum enclosing circle O(n)
        Circle mec() const {
            return welzl(ps);
        }

        friend istream& operator>>(istream& in, Polygon& q) {
            for (P& p: q.ps) in >> p;
            return in;
        }

        friend ostream& operator<<(ostream& out, const Polygon& q) {
            if (q.ps.empty()) return out;
            for (const P& p: q.ps) out << p << ", ";
            out << q.ps[0];
            return out;
        }
    };

    // Monotone chain, sorted by lexi (x, y), from CP
    Polygon convexHullMonotone(vec<P> a, bool include_collinear = false) {
        if (a.size() == 1) return Polygon{a};

        sort(all(a));
        P p1 = a[0], p2 = a.back();
        vector<P> up, down;
        up.push_back(p1);
        down.push_back(p1);

        for (int i = 1; i < (int)a.size(); i++) {
            if (i == a.size() - 1 or cw(p1, a[i], p2, include_collinear)) {
                while (up.size() >= 2 and not cw(up[up.size()-2], up[up.size()-1], a[i], include_collinear))
                    up.pop_back();
                up.push_back(a[i]);
            }
            if (i == a.size() - 1 or ccw(p1, a[i], p2, include_collinear)) {
                while (down.size() >= 2 and not ccw(down[down.size()-2], down[down.size()-1], a[i], include_collinear))
                    down.pop_back();
                down.push_back(a[i]);
            }
        }

        if (include_collinear and up.size() == a.size()) {
            reverse(all(a));
            return Polygon{a};
        }
        a.clear();
        for (const auto & i : up) a.push_back(i);
        for (int i = (int)down.size() - 2; i > 0; i--) a.push_back(down[i]);
        return Polygon{a};
    }

    // Based on CSES
    Polygon convexHullMonotoneSm(vec<P> p) {
        sort(all(p));
        int n = (int)p.size();
        vec<P> res;
        int s = 0;

        for (int t = 0; t < 2; t++) {
            for (int i = 0; i < n; i++) {
                while (res.size() - s >= 2) {
                    P p1 = res[res.size() - 2];
                    P p2 = res[res.size() - 1];
                    if (p1.cross(p2, p[i]) <= 0) break;
                    res.pop_back();
                }

                res.push_back(p[i]);
            }
            res.pop_back();
            s = (int)res.size();
            reverse(all(p));
        }

        return Polygon{ res };
    }

    // cw order starting from 12 o'clock, if on same hour, shortest distance from center.
    void radialSort(vec<P>& v) {
        P c;
        for (P& p: v) c += p;
        c /= v.size();
        auto less = [&](const P& a, const P& b) -> bool {
            if (lse(0, a.x - c.x) and ls(b.x - c.x, 0)) return true;
            if (ls(a.x - c.x, 0) and lse(0, b.x - c.x)) return false;
            if (eq(a.x - c.x, 0) and eq(b.x - c.x, 0)) {
                if (lse(0, a.y - c.y) or lse(0, b.y - c.y)) return ls(b.y, a.y);
                return ls(a.y, b.y);
            }
            ftype det = (a - c).cross(b - c);
            if (ls(det, 0)) return true;
            if (ls(0, det)) return false;

            // points a and b are on the same line from the center
            // check which point is closer to the center
            ftype d1 = (a.x - c.x) * (a.x - c.x) + (a.y - c.y) * (a.y - c.y),
                  d2 = (b.x - c.x) * (b.x - c.x) + (b.y - c.y) * (b.y - c.y);

            return ls(d2, d1);
        };
        sort(all(v), less);
    }

    ftype ccw(const P& a, const P& b, const P& ref = {0, 0}) {
        return (a - ref).cross(b - ref);
    }

    struct AngleCompare {
        const P origin;
        explicit AngleCompare(const P& origin = {}): origin{origin} {}
        bool operator()(const P& l, const P& r) const {
            ftype tmp = ccw(l, r, origin);
            return eq(tmp, 0) ? ls(l.dist2(origin), r.dist2(origin)) : ls(tmp, 0);
        }
    };

    // cw order starting from 12 o'clock, if on same hour, one closest to reference is chosen
    // works with integers.
    template<class I>
    void sortByAngle(I first, I last, const P& ref = {0, 0}) {
        first = partition(first, last, [&ref](const decltype(*first)& p) {
            return p == ref;
        });
        auto pivot = partition(first, last, [&ref](const decltype(*first)& p) {
            return ref < p;
        });
        AngleCompare acmp(ref);
        sort(first, pivot, acmp);
        sort(pivot, last, acmp);
    }

    int maxByManhattan(const vec2d<int>& p) {
        // O(n * d * 2 ^ n)
        int res = 0;
        const int static d = 2;
        for (int msk = 0; msk < (1ll << d); msk++) {
            int mx = LLONG_MIN, mn = LLONG_MAX;
            for (const auto& r : p) {
                int cur = 0;
                for (int j = 0; j < d; j++) {
                    if (msk & (1ll << j)) cur += r[j];
                    else cur -= r[j];
                }
            }
            res = max(res, mx - mn);
        }
        return res;
    }

    // Sorts points from left to right in the given direction
    // (-1, 0) is equivalent to (1, 0)
    void sortInDir(vec<P>& v, const V& d) {
        sort(all(v), [&](const P& a, const P& b) -> bool {
            const P& ap = d.project(a), bp = d.project(b);
            return ap == bp ? a < b : ap < bp;
        });
    }

    vec<pii> fortunes(const vec<P>& p) {
        using namespace vornoi;
        auto f = fortune(p);
        f.solve();
        return f.edges;
    }

    ostream& operator<<(ostream& out, const vec<Triangle>& a) {
        for (auto& t: a) {
            out << get<0>(t) << ", " << get<1>(t) << ", " << get<2>(t) << ", " << get<0>(t) << nl;
        }
        return out;
    }

    // id of the vertex having maximum dot product with z
    // polygon must need to be convex
    // top - upper right vertex
    // for minimum dot product negate z and return -dot(z, p[id])
    int extreme_vertex(vector<P>& p, const P& z, const int top) { // O(log n)
        int n = p.size();
        if (n == 1) return 0;
        ftype ans = p[0].dot(z); int id = 0;
        if (ls(ans, p[top].dot(z))) ans = p[top].dot(z), id = top;
        int l = 1, r = top - 1;
        while (l < r) {
            int mid = (l + r) >> 1ll;
            if (lse(p[mid].dot(z), p[mid + 1].dot(z))) l = mid + 1;
            else r = mid;
        }
        if (ls(ans, p[l].dot(z))) ans = p[l].dot(z), id = l;
        l = top + 1, r = n - 1;
        while (l < r) {
            int mid = (l + r) >> 1ll;
            if (lse(p[mid].dot(z), p[(mid + 1) % n].dot(z))) l = mid + 1;
            else r = mid;
        }
        l %= n;
        if (ls(ans, p[l].dot(z))) ans = p[l].dot(z), id = l;
        return id;
    }

    // given a convex polygon p, and a line ab and the top vertex of the polygon
    // returns the intersection of the line with the polygon
    // it returns the indices of the edges of the polygon that are intersected by the line
    // so if it returns i, then the line intersects the edge (p[i], p[(i + 1) % n])
    array<int, 2> convex_line_intersection(vector<P>& p, P a, P b, int top) {
        int end_a = extreme_vertex(p, (a - b).normal_r(), top);
        int end_b = extreme_vertex(p, (b - a).normal_r(), top);
        auto cmp_l = [&](int i) { return orientation(a, p[i], b); };
        if (cmp_l(end_a) < 0 || cmp_l(end_b) > 0)
            return {-1, -1}; // no intersection
        array<int, 2> res;
        for (int i = 0; i < 2; i++) {
            int lo = end_b, hi = end_a, n = p.size();
            while ((lo + 1) % n != hi) {
                int m = ((lo + hi + (lo < hi ? 0 : n)) / 2) % n;
                (cmp_l(m) == cmp_l(end_b) ? lo : hi) = m;
            }
            res[i] = (lo + !cmp_l(hi)) % n;
            swap(end_a, end_b);
        }
        if (res[0] == res[1]) return {res[0], -1}; // touches the vertex res[0]
        if (!cmp_l(res[0]) && !cmp_l(res[1]))
            switch ((res[0] - res[1] + (int)p.size() + 1) % p.size()) {
                case 0: return {res[0], res[0]}; // touches the edge (res[0], res[0] + 1)
                case 2: return {res[1], res[1]}; // touches the edge (res[1], res[1] + 1)
            }
        return res; // intersects the edges (res[0], res[0] + 1) and (res[1], res[1] + 1)
    }

    pair<P, int> point_poly_tangent(vector<P>& p, P Q, int dir, int l, int r) {
        while (r - l > 1) {
            int mid = (l + r) >> 1;
            bool pvs = orientation(Q, p[mid], p[mid - 1]) != -dir;
            bool nxt = orientation(Q, p[mid], p[mid + 1]) != -dir;
            if (pvs and nxt) return {p[mid], mid};
            if (not (pvs or nxt)) {
                const auto& p1 = point_poly_tangent(p, Q, dir, mid + 1, r);
                const auto& p2 = point_poly_tangent(p, Q, dir, l, mid - 1);
                return orientation(Q, p1.first, p2.first) == dir ? p1 : p2;
            }
            if (!pvs) {
                if (orientation(Q, p[mid], p[l]) == dir)  r = mid - 1;
                else if (orientation(Q, p[l], p[r]) == dir) r = mid - 1;
                else l = mid + 1;
            }
            if (!nxt) {
                if (orientation(Q, p[mid], p[l]) == dir)  l = mid + 1;
                else if (orientation(Q, p[l], p[r]) == dir) r = mid - 1;
                else l = mid + 1;
            }
        }
        pair<P, int> ret = {p[l], l};
        for (int i = l + 1; i <= r; i++) ret = orientation(Q, ret.first, p[i]) != dir ? make_pair(p[i], i) : ret;
        return ret;
    }
    // (ccw, cw) tangents from a point that is outside this convex polygon
    // returns indexes of the points
    // ccw means the tangent from Q to that point is in the same direction as the polygon ccw direction
    pair<int, int> tangents_from_point_to_polygon(vector<P> &p, P Q) {
        int ccw = point_poly_tangent(p, Q, 1, 0, (int)p.size() - 1).second;
        int cw = point_poly_tangent(p, Q, -1, 0, (int)p.size() - 1).second;
        return make_pair(ccw, cw);
    }

    // minimum distance from point c to segment ab
    double dist_from_point_to_seg(const P& a, const P& b, const P& c) {
        return V(a, b).dist(c);
    }

    // minimum distance from point c to line through a and b
    double dist_from_point_to_line(const P& a, const P& b, const P& c) {
        return fabs((b - a).cross(c - a) / (b - a).norm());
    }

    // minimum distance from a point to a convex polygon
    // it assumes point lie strictly outside the polygon
    double dist_from_point_to_polygon(vector<P> &p, P z) {
        double ans = inf;
        int n = p.size();
        if (n <= 3) {
            for(int i = 0; i < n; i++) ans = min(ans, dist_from_point_to_seg(p[i], p[(i + 1) % n], z));
            return ans;
        }
        auto [r, l] = tangents_from_point_to_polygon(p, z);
        if(l > r) r += n;
        while (l < r) {
            int mid = (l + r) >> 1;
            double left = p[mid % n].dist2(z), right = p[(mid + 1) % n].dist2(z);
            ans = std::min({ans, left, right});
            if(left < right) r = mid;
            else l = mid + 1;
        }
        ans = sqrt(ans);
        ans = min(ans, dist_from_point_to_seg(p[l % n], p[(l + 1) % n], z));
        ans = min(ans, dist_from_point_to_seg(p[l % n], p[(l - 1 + n) % n], z));
        return ans;
    }
    // minimum distance from convex polygon p to line ab
    // returns 0 is it intersects with the polygon
    // top - upper right vertex
    double dist_from_polygon_to_line(vector<P> &p, P a, P b, int top) { //O(log n)
        P orth = (b - a).normal_r();
        if (orientation(a, b, p[0]) > 0) orth = (a - b).normal_r();
        int id = extreme_vertex(p, orth, top);
        if (ls(0, (p[id] - a).dot(orth))) return 0.0; //if orth and a are in the same half of the line, then poly and line intersects
        return dist_from_point_to_line(a, b, p[id]); //does not intersect
    }
    // minimum distance from a convex polygon to another convex polygon
    // the polygon doesnot overlap or touch
    // tested in https://toph.co/p/the-wall
    double dist_from_polygon_to_polygon(vector<P> &p1, vector<P> &p2) { // O(n log n)
        double ans = inf;
        for (int i = 0; i < p1.size(); i++) {
            ans = min(ans, dist_from_point_to_polygon(p2, p1[i]));
        }
        for (int i = 0; i < p2.size(); i++) {
            ans = min(ans, dist_from_point_to_polygon(p1, p2[i]));
        }
        return ans;
    }
    // maximum distance from a convex polygon to another convex polygon
    double maximum_dist_from_polygon_to_polygon(vector<P> &u, vector<P> &v){ //O(n)
        int n = (int)u.size(), m = (int)v.size();
        double ans = 0;
        if (n < 3 || m < 3) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) ans = max(ans, u[i].dist2(v[j]));
            }
            return sqrt(ans);
        }
        if (u[0].x > v[0].x) swap(n, m), swap(u, v);
        int i = 0, j = 0, step = n + m + 10;
        while (j + 1 < m && v[j].x < v[j + 1].x) j++ ;
        while (step--) {
            if ((u[(i + 1)%n] - u[i]).cross(v[(j + 1)%m] - v[j]) >= 0) j = (j + 1) % m;
            else i = (i + 1) % n;
            ans = max(ans, u[i].dist2(v[j]));
        }
        return sqrt(ans);
    }

    // calculates the area of the union of n polygons (not necessarily convex).
    // the points within each polygon must be given in CCW order.
    // complexity: O(N^2), where N is the total number of points
    double rat(const P& a, const P& b, const P& p) {
        return !sign(a.x - b.x) ? (p.y - a.y) / (b.y - a.y) : (p.x - a.x) / (b.x - a.x);
    };
    double polygon_union(vector<vector<P>> &p) {
        int n = p.size();
        double ans=0;
        for(int i = 0; i < n; ++i) {
            for (int v = 0; v < (int)p[i].size(); ++v) {
                P a = p[i][v], b = p[i][(v + 1) % p[i].size()];
                vector<pair<double, int>> segs;
                segs.emplace_back(0,  0), segs.emplace_back(1,  0);
                for(int j = 0; j < n; ++j) {
                    if(i != j) {
                        for(size_t u = 0; u < p[j].size(); ++u) {
                            P c = p[j][u], d = p[j][(u + 1) % p[j].size()];
                            int sc = sign((b - a).cross(c - a)), sd = sign((b - a).cross(d - a));
                            if(!sc && !sd) {
                                if(sign((b - a).dot(d - c)) > 0 && i > j) {
                                    segs.emplace_back(rat(a, b, c), 1), segs.emplace_back(rat(a, b, d),  -1);
                                }
                            }
                            else {
                                double sa = (d - c).cross(a - c), sb = (d - c).cross(b - c);
                                if(sc >= 0 && sd < 0) segs.emplace_back(sa / (sa - sb), 1);
                                else if(sc < 0 && sd >= 0) segs.emplace_back(sa / (sa - sb),  -1);
                            }
                        }
                    }
                }
                sort(segs.begin(),  segs.end());
                double pre = min(max(segs[0].first, 0.0), 1.0), now, sum = 0;
                int cnt = segs[0].second;
                for(int j = 1; j < segs.size(); ++j) {
                    now = min(max(segs[j].first, 0.0), 1.0);
                    if (!cnt) sum += now - pre;
                    cnt += segs[j].second;
                    pre = now;
                }
                ans += a.cross(b) * sum;
            }
        }
        return ans * 0.5;
    }

    // returns the area of the intersection of the circle with center c and radius r
    // and the triangle formed by the points c, a, b
    double _triangle_circle_intersection(P c, double r, P a, P b) {
        double sd1 = c.dist2(a), sd2 = c.dist2(b);
        if(sd1 > sd2) swap(a, b), swap(sd1, sd2);
        double sd = a.dist2(b);
        double d1 = sqrtl(sd1), d2 = sqrtl(sd2), d = sqrt(sd);
        double x = abs(sd2 - sd - sd1) / (2 * d);
        double h = sqrtl(sd1 - x * x);
        if(r >= d2) return h * d / 2;
        double area = 0;
        if(sd + sd1 < sd2) {
            if(r < d1) area = r * r * (acos(h / d2) - acos(h / d1)) / 2;
            else {
                area = r * r * ( acos(h / d2) - acos(h / r)) / 2;
                double y = sqrtl(r * r - h * h);
                area += h * (y - x) / 2;
            }
        }
        else {
            if(r < h) area = r * r * (acos(h / d2) + acos(h / d1)) / 2;
            else {
                area += r * r * (acos(h / d2) - acos(h / r)) / 2;
                double y = sqrtl(r * r - h * h);
                area += h * y / 2;
                if(r < d1) {
                    area += r * r * (acos(h / d1) - acos(h / r)) / 2;
                    area += h * y / 2;
                }
                else area += h * x / 2;
            }
        }
        return area;
    }
    // intersection between a simple polygon and a circle
    double polygon_circle_intersection(vector<P> &v, P p, double r) {
        int n = v.size();
        double ans = 0.00;
        P org = {0, 0};
        for(int i = 0; i < n; i++) {
            int x = orientation(p, v[i], v[(i + 1) % n]);
            if(x == 0) continue;
            double area = _triangle_circle_intersection(org, r, v[i] - p, v[(i + 1) % n] - p);
            if (x < 0) ans -= area;
            else ans += area;
        }
        return abs(ans);
    }
}
