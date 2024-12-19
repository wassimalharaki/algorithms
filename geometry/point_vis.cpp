#include "./polygon.cpp"
using Poly = vec<V>;
enum Orientation {
    left_turn = 1,
    right_turn = -1,
    collinear = 0
};
int getOrientation(const P& a, const P& b, const P& c) {
    return sign((b - a).cross(c - a));
}
struct Ray {
    P org, dir;
    Ray(): org{}, dir{} {};
    Ray(const P& orig, const P& dir): org{orig}, dir{dir} {}
    bool intersects(const V& seg, P& outp) const {
        auto ao = org - seg.s;
        auto ab = seg.d - seg.s;
        auto det = ab.cross(dir);
        if (eq(det, 0)) {
            auto abo = getOrientation(seg.s, seg.d, org);
            if (abo != Orientation::collinear) return false;
            auto dist_a = ao.dot(dir);
            auto dist_b = (org - seg.d).dot(dir);

            if (ls(0, dist_a) and ls(0, dist_b)) return false;
            else if (ls(dist_a, 0) != ls(0, dist_b)) outp = org;
            else if (ls(dist_b, dist_a)) outp = seg.s;
            else outp = seg.d;
            return true;
        }
        auto u = ao.cross(dir) / det;
        if (ls(u, 0) or ls(1, u)) return false;
        auto t = -ab.cross(ao) / det;
        outp = org + t * dir;
        return lse(0, t);
    }
};
struct LineSegmentDistComparer {
    using segment_type = V;
    P origin;
    explicit LineSegmentDistComparer(const P& origin): origin(origin) {}

    /** Check whether the line segment x is closer to the origin than the
     * line segment y.
     * @param x line segment: left hand side of the comparison operator
     * @param y line segment: right hand side of the comparison operator
     * @return true iff x < y (x is closer than y)
     */
    bool operator()(const segment_type& x, const segment_type& y) const {
        auto a = x.s, b = x.d;
        auto c = y.s, d = y.d;

        assert(
                getOrientation(origin, a, b) != Orientation::collinear &&
                "AB must not be collinear with the origin.");
        assert(
                getOrientation(origin, c, d) != Orientation::collinear &&
                "CD must not be collinear with the origin.");

        // sort the endpoints so that if there are common endpoints,
        // it will be a and c
        if (b == c or b == d)
            std::swap(a, b);
        if (a == d)
            std::swap(c, d);

        // cases with common endpoints
        if (a == c)
        {
            auto oad = getOrientation(origin, a, d);
            auto oab = getOrientation(origin, a, b);
            if (b == d || oad != oab)
                return false;
            return getOrientation(a, b, d) != getOrientation(a, b, origin);
        }

        // cases without common endpoints
        auto cda = getOrientation(c, d, a);
        auto cdb = getOrientation(c, d, b);
        if (cdb == Orientation::collinear && cda == Orientation::collinear) {
            return (origin - a).norm() < (origin - c).norm();
        }
        else if (cda == cdb ||
                 cda == Orientation::collinear ||
                 cdb == Orientation::collinear)
        {
            auto cdo = getOrientation(c, d, origin);
            return cdo == cda || cdo == cdb;
        }
        else
        {
            auto abo = getOrientation(a, b, origin);
            return abo != getOrientation(a, b, c);
        }
    }
};
struct AngleComparer {
    P vertex;
    explicit AngleComparer(const P& origin): vertex(origin) {}

    bool operator()(const P& a, const P& b) const {
        auto is_a_left = ls(a.x, vertex.x);
        auto is_b_left = ls(b.x, vertex.x);
        if (is_a_left != is_b_left)
            return is_b_left;

        if (eq(a.x, vertex.x) && eq(b.x, vertex.x)) {
            if (!ls(a.y, vertex.y) ||
                !ls(b.y, vertex.y))
            {
                return ls(b.y, a.y);
            }
            return ls(a.y, b.y);
        }

        auto oa = a - vertex;
        auto ob = b - vertex;
        auto det = oa.cross(ob);
        if (eq(det, 0.f)) {
            return ls(oa.norm(), ob.norm());
        }
        return ls(det, 0);
    }
};
struct VisibilityEvent {
    enum EventType {
        start_vertex, end_vertex
    };
    EventType type{};
    V segment;
    VisibilityEvent() = default;
    VisibilityEvent(EventType type, V segment): type(type), segment(std::move(segment)) {}
    [[nodiscard]] const auto& point() const { return segment.s; }
};
Polygon visibility_polygon(const P& point, const Poly& p) {
    using segment_type = V;
    using event_type = VisibilityEvent;
    using segment_comparer_type = LineSegmentDistComparer;
    segment_comparer_type cmp_dist{ point };
    std::set<segment_type, segment_comparer_type> state{ cmp_dist };
    std::vector<event_type> events;
    for (auto& segment: p) {
        // Sort line segment endpoints and add them as events
        // Skip line segments collinear with the point
        auto pab = getOrientation(point, segment.s, segment.d);
        if (pab == Orientation::collinear)
        {
            continue;
        }
        else if (pab == Orientation::right_turn)
        {
            events.emplace_back(event_type::start_vertex, segment);
            events.emplace_back(
                    event_type::end_vertex,
                    segment_type{ segment.d, segment.s });
        }
        else
        {
            events.emplace_back(
                    event_type::start_vertex,
                    segment_type{ segment.d, segment.s });
            events.emplace_back(event_type::end_vertex, segment);
        }

        // Initialize state by adding line segments that are intersected
        // by vertical ray from the point
        auto a = segment.s, b = segment.d;
        if (a.x > b.x)
            std::swap(a, b);

        auto abp = getOrientation(a, b, point);
        if (abp == Orientation::right_turn &&
            (eq(b.x, point.x) ||
             (a.x < point.x && point.x < b.x)))
        {
            state.insert(segment);
        }
    }

    // sort events by angle
    AngleComparer cmp_angle{ point };
    std::sort(events.begin(), events.end(), [&cmp_angle](auto&& a, auto&& b)
    {
        // if the points are eq, sort end vertices first
        if (a.point() == b.point())
            return a.type == event_type::end_vertex &&
                   b.type == event_type::start_vertex;
        return cmp_angle(a.point(), b.point());
    });

    // find the visibility polygon
    vec<P> vertices;
    for (auto&& event : events)
    {
        if (event.type == event_type::end_vertex)
            state.erase(event.segment);

        if (state.empty())
        {
            vertices.push_back(event.point());
        }
        else if (cmp_dist(event.segment, *state.begin()))
        {
            // Nearest line segment has changed
            // Compute the intersection point with this segment
            P intersection;
            Ray ray{ point, event.point() - point };
            auto nearest_segment = *state.begin();
            auto intersects = ray.intersects(nearest_segment, intersection);
            assert(intersects &&
                           "Ray intersects line segment L iff L is in the state");

            if (event.type == event_type::start_vertex)
            {
                vertices.push_back(intersection);
                vertices.push_back(event.point());
            }
            else
            {
                vertices.push_back(event.point());
                vertices.push_back(intersection);
            }
        }

        if (event.type == event_type::start_vertex)
            state.insert(event.segment);
    }

    // remove collinear points
    auto top = vertices.begin();
    for (auto it = vertices.begin(); it != vertices.end(); ++it)
    {
        auto prev = top == vertices.begin() ? vertices.end() - 1 : top - 1;
        auto next = it + 1 == vertices.end() ? vertices.begin() : it + 1;
        if (getOrientation(*prev, *it, *next) != Orientation::collinear)
            *top++ = *it;
    }
    vertices.erase(top, vertices.end());
    for (auto& s: p) {
        if (s.s == point or s.d == point) {
            vertices.push_back(point);
            break;
        }
    }
    return Polygon{vertices};
}
