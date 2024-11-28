#include <bits/stdc++.h>
using namespace std;

namespace geometry {

    /// Underlying type for all computations
    using float_t = long double;

    /// Permissible error; in other words, limiting accuracy
    inline constexpr float_t EPS = 1e-9l;

    /// \cond
    inline constexpr float_t EPS_INC = 100 * EPS;
    /// \endcond

} // namespace geometry

namespace geometry {

    /**
     * @brief Simple point struct to store the x and y coordinates of point entities
     */
    struct point_t {
        /// The x coordinate of the point
        float_t x;

        /// The y coordinate of the point
        float_t y;

        /**
         * @brief Overloading the == operator for point_t
         *
         * All floating point comparisons are done within a neighbourhood of `geometry::EPS`.
         *
         * @param other The other point this is being compared to
         * @return `true` if both the x and y coordinates match with max permissible error `geometry::EPS`
         * @return `false` otherwise
         */
        bool operator == (const point_t &other) const;
    };

} // namespace geometry

namespace geometry {

    /**
     * @brief Segment struct which encapsulates information about a segment
     * @warning \f$ p \le q \f$ must hold for the algorithm to work (comparing on the tuple `<x, y>`)
     */
    struct segment_t {
        /// The point where the segment begins
        point_t p;

        /// The point where the segment ends
        point_t q;

        /// The index of the segment in input
        size_t seg_id;

        /**
         * @brief Evaluates the y coordinate of a point on the segment given its x coordinate
         *
         * @param x The x coordinate of the point to be found on the segment
         * @return `y` The corresponding y coordinate of the desired point
         */
        float_t eval_y(float_t x) const;

        /**
         * @brief  Overloading the < operator overload for segment_t
         *
         * Required by (`std::less<T>`, the defalt comparator of) the segment ordering BBST
         * to make comparisons between events and establish an ordering.
         *
         * Compares y coordinates of segments by calling `segment_t::eval_y()` with the current x coordinate of the sweepline. <br>
         * All floating point comparisons are done within a neighbourhood of `geometry::EPS`.
         *
         * @param s The other segment to compare to
         * @return `true` if the y coordinate corresponding to `sweepline::sweeplineX` is lesser than the other segment's
         * @return `false` otherwise
         */
        bool operator < (const segment_t &s) const;
    };

    /**
     * @brief Checks if two segments intersect in one dimension
     *
     * @param l1 `p.x` (or `p.y`) of the first segment
     * @param r1 `q.x` (or `q.y`) of the first segment
     * @param l2 `p.x` (or `p.y`) of the second segment
     * @param r2 `q.x` (or `q.y`) of the second segment
     * @return `true` if the segments intersect in one dimension
     * @return `false` otherwise
     */
    bool can_intersect_1d(float_t l1, float_t r1, float_t l2, float_t r2);

    /**
     * @brief Computes the sign of the cross product of three points
     *
     * @param a The first point
     * @param b The second point
     * @param c The third point
     * @return `-1` if the cross product is negative
     * @return `0`  if the cross product is zero
     * @return `+1` if the cross product is positive
     */
    int cross_prod(const point_t &a, const point_t &b, const point_t &c);

    /**
     * @brief Checks if two segments intersect
     *
     * Necessary and sufficient conditions for intersection:
     * 1. Segments must intersect one dimensionally on both x and y components.
     * 2. Cross product of the end points of each segment with each end point of
     * the other segment must not have the same sign.
     *
     * @param a The first segment
     * @param b The second segment
     * @return `true` if the segments intersect
     * @return `false` otherwise
     */
    bool is_intersecting(const segment_t &a, const segment_t &b);

    /**
     * @brief Computes the point of intersection of two segments
     * @pre The segments must intersect.
     *
     * @param a The first segment
     * @param b The second segment
     * @return The point of intersection of the two segments
     */
    point_t intersection_point(const segment_t &a, const segment_t &b);

} // namespace geometry

namespace sweepline {

    /**
     * @brief Event struct which encapsulates information about an event
     *
     * @note Multiple events having the same event point may exist simultaneously
     * within the event queue, one for each segment containing its `seg_id`.
     *
     */
    struct event_t {
        /// The point where the event occurs
        geometry::point_t p;

        /// An enum for the type of event
        enum type {
            begin,    ///< denotes when a segment starts
            interior, ///< denotes an interior point of a segment
            end       ///< denotes when a segment ends
        } tp;   ///< The type of the event

        /// The id of the segment
        size_t seg_id;

        /**
         * @brief Default constructor
         */
        event_t() = default;

        /**
         * @brief Parameterized constructor
         *
         * @param p The point where the event occurs
         * @param tp The type of the event
         * @param seg_id The id of the segment
         */
        event_t(const geometry::point_t &p, type tp, size_t seg_id)
                : p(p), tp(tp), seg_id(seg_id) {}

        /**
         * @brief Overloading the < operator for event_t
         *
         * Required by (`std::less<T>`, the defalt comparator of) the event queue BBST
         * to make comparisons between events and establish an ordering.
         *
         * Compares on, in decreasing priority, the tuple `<p.x, p.y, p.seg_id>`. <br>
         * All floating point comparisons are done within a neighbourhood of `geometry::EPS`.
         *
         * @param e The other event to be compared to
         * @return `true` if it compares less than the other event
         * @return `false` otherwise
         */
        bool operator < (const event_t &e) const;

    };

} // namespace sweepline

namespace sweepline {

    /// A global variable which stores the current X coordinate of the vertical sweepline
    extern geometry::float_t sweeplineX;

    /**
     * @brief Simple struct to bind together information on an intersection of two or more segments
     */
    struct intersection_t {
        geometry::point_t pt;   ///< The point of intersection of some segments
        std::vector<size_t> segments;   ///< A list of segments (their ids) which intersect at this point
    };

    /**
     * @brief Finds which segments intersect at which points and returns all such intersections
     * @pre \f$ p \le q \f$ must hold for each line segment in the list.
     * @pre No two line segments should coincide with each other either fully or partially.
     *
     * Takes a list of line segments whose end points are assumed to be in sorted order (left to right).
     *
     * Uses the Bentley Ottman algorithm to compute all line segment intersections efficiently in \f$ \mathcal{O}((n+k)logn) \f$, <br>
     * where \f$ n \f$ is the number of line segments and \f$ k \f$ is the number of intersection points.
     *
     * Returns a list of intersections, i.e. pairs of points and the corresponding
     * indices (1-based) of segments which intersect at that point.
     *
     * Calls `solver::solve()` and returns the result.
     *
     * @param line_segments The list of input line segments
     * @param verbose The `utils::args::verbose` flag
     * @param enable_color The `utils::args::enable_color` flag
     * @return `std::vector<intersection_t>` A list of all intersections
     */
    std::vector<intersection_t> find_intersections(
            const std::vector<geometry::segment_t> &line_segments,
            bool verbose = false,
            bool enable_color = true
    );

    template <typename T>
//    using bbst = BBST::red_black_tree<T>;   ///< Type alias for the underlying BBST used. Works with std::set in exactly the same way as well.
     using bbst = std::set<T>; // works with std::set in exactly the same way (don't forget to #include <set>)

    /**
     * @brief A utility class instantiated by `find_intersections()`
     * to manage the data structures and implement the algorithm
     * to find all intersection points.
     *
     */
    class solver {
        bool verbose;                                     ///< The `utils::args::verbose` flag
        std::vector<geometry::segment_t> line_segments;   ///< The list of input line segments
        std::vector<sweepline::intersection_t> result;    ///< The list of intersections that will be returned

        bbst<sweepline::event_t> event_queue;             ///< The event queue, implemented as a BBST of events
        bbst<geometry::segment_t> seg_ordering;           ///< The status queue, or segment ordering, implemented as a BBST of segments
        std::vector<geometry::segment_t> vertical_segs;   ///< A list of line segments with slope parallel to the sweepline (vertical) that will be handled separately

    public:

        /**
         * @brief Constructor
         *
         * @param line_segments The list of input line segments
         * @param verbose The `utils::args::verbose` flag
         * @param enable_color The `utils::args::enable_color` flag
         */
        solver(const std::vector<geometry::segment_t> &line_segments, bool verbose, bool enable_color);

        /**
         * @brief Finds which segments intersect at which points and returns all such intersections
         *
         * @pre \f$ p \le q \f$ must hold for each line segment in the list.
         * @pre No two line segments should coincide with each other either fully or partially.
         *
         * Takes a list of line segments whose end points are assumed to be in sorted order (left to right).
         *
         * Uses the Bentley Ottman algorithm to compute all line segment intersections efficiently in \f$ \mathcal{O}((n+k)logn) \f$, <br>
         * where \f$ n \f$ is the number of line segments and \f$ k \f$ is the number of intersection points.
         *
         * Returns a list of intersections, i.e. pairs of points and the corresponding
         * indices (1-based) of segments which intersect at that point.
         *
         * @return `std::vector<intersection_t>` A list of all intersections
         */
        std::vector<sweepline::intersection_t> solve();

    private:
        // Implementation

        /**
         * @brief Initializes the `solver::event_queue` by inserting the end points of the `solver::line_segments` and populates `solver::vertical_segs` with vertical segments
         * @pre \f$ p \le q \f$ must hold for each line segment in the list.
         *
         * Iterates over the segments in `solver::line_segments`
         * 1. Vertical segments are added to `solver::vertical_segs`
         * 2. For all other segments, the begin and end points are inserted into the `solver::event_queue` as begin and end events respectively
         *
         */
        void init_event_queue();

        /**
         * @brief Finds intersections between pairs of vertical line segments
         *
         * There can be at most two non-coinciding and intersecting vertical line segments.
         * Such intersections are found by checking if adjacent pairs of segments in sorted order (by x and then y) intersect.
         *
         */
        void find_vertical_vertical_intersections();

        /**
         * @brief Finds intersections between a vertical line segment and one or more non-vertical segments
         *
         * Iterates over vertical line segments (ordered by x) whose `x` coordinate lies before the new sweepline.
         * Non-vertical segments which intersect with a particular vertical segment are found by querying `solver::seg_ordering`.
         *
         * @param max_vsegx The `x` coordinate of the next event, i.e. the new sweepline.
         */
        void find_vertical_nonvertical_intersections(geometry::float_t max_vsegx);

        /**
         * @brief Gets the active segments with an event at the point currently being processed
         *
         * @param top One of the events at the point currently being processed
         * @return `std::array<std::vector<size_t>, 3>` Three arrays of active segment indices corresponding to `event_t::type`
         */
        std::array<std::vector<size_t>, 3> get_active_segs(sweepline::event_t top);

        /**
         * @brief Updates the status queue `solver::seg_ordering` after processing the event point
         *
         * * Removes all segments with `event_t::type::end` and `event_t::type::interior` events
         * * Increments the sweepline by a very small amount, just past the intersection point
         * * Inserts all segments with `event_t::type::begin` events and reinserts the removed
         * segments with `event_t::type::interior` events in the correct order
         *
         * @param active_segs The active segments with an event at the point currently being processed
         */
        void update_segment_ordering(const std::array<std::vector<size_t>, 3> &active_segs);

        /**
         * @brief Tests for new event points after updating `solver::seg_ordering` in the case when no new segments are inserted
         *
         * If no segments were newly inserted, the immediate left and right neighbours
         * of the deleted set of segments become adjacent candidates for intersection.
         *
         * @param cur The current point being processed
         */
        void handle_no_newly_inserted(geometry::point_t cur);

        /**
         * @brief Tests for new event points after updating `solver::seg_ordering` in the case when some new segments are inserted
         *
         * If some segments were newly inserted,
         * the left and right extremes among the set of newly inserted segments
         * must be checked for intersection with their immediate left and right neighbours respectively.
         *
         */
        void handle_extremes_of_newly_inserted();

        /**
         * @brief Reports an intersection between teo or more (non-vertical) line segments
         *
         * @param cur The point of intersection
         * @param active_segs The line segments that intersect at \a cur
         */
        void report_intersection(geometry::point_t cur, std::array<std::vector<size_t>, 3> &&active_segs);

        /**
         * @brief Merge intersections which have the same point
         *
         */
        void merge_intersection_points();

        /// \cond
        size_t vert_idx = 0;
        geometry::float_t max_y, min_y;
        /// \endcond
    };

} // namespace sweepline

bool geometry::point_t::operator == (const point_t &other) const {
    return std::fabs(x - other.x) < EPS and std::fabs(y - other.y) < EPS;
}

namespace sweepline {
    /// A global variable which stores the current X coordinate of the vertical sweepline
    extern geometry::float_t sweeplineX;
}

geometry::float_t geometry::segment_t::eval_y(geometry::float_t x) const {
    return std::fabs(p.x - q.x) < EPS? p.y
                                     : p.y + (q.y - p.y) * (x - p.x) / (q.x - p.x);
}

bool geometry::segment_t::operator < (const segment_t &other) const {
    return eval_y(sweepline::sweeplineX) < other.eval_y(sweepline::sweeplineX) - EPS;
}

bool geometry::can_intersect_1d(
        geometry::float_t l1, geometry::float_t r1, geometry::float_t l2, geometry::float_t r2
) {
    if(l1 > r1) std::swap(l1, r1);
    if(l2 > r2) std::swap(l2, r2);
    return std::max(l1, l2) <= std::min(r1, r2) + EPS;
}

int geometry::cross_prod(
        const geometry::point_t &a, const geometry::point_t &b, const geometry::point_t &c
) {
    geometry::float_t s = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    return std::fabs(s) < EPS ? 0 : s > 0 ? +1 : -1;
}

bool geometry::is_intersecting(const segment_t &a, const segment_t &b) {
    return can_intersect_1d(a.p.x, a.q.x, b.p.x, b.q.x)
           and can_intersect_1d(a.p.y, a.q.y, b.p.y, b.q.y)
           and cross_prod(a.p, a.q, b.p) * cross_prod(a.p, a.q, b.q) <= 0
           and cross_prod(b.p, b.q, a.p) * cross_prod(b.p, b.q, a.q) <= 0;
}

geometry::point_t geometry::intersection_point(const segment_t &a, const segment_t &b) {
    geometry::float_t A1 = a.p.y - a.q.y;
    geometry::float_t B1 = a.q.x - a.p.x;
    geometry::float_t C1 = (a.p.y - a.q.y) * (-a.p.x) + (a.q.x - a.p.x) * (-a.p.y);

    geometry::float_t A2 = b.p.y - b.q.y;
    geometry::float_t B2 = b.q.x - b.p.x;
    geometry::float_t C2 = (b.p.y - b.q.y) * (-b.p.x) + (b.q.x - b.p.x) * (-b.p.y);

    return geometry::point_t {(B1 * C2 - B2 * C1) / (A1 * B2 - A2 * B1), (C1 * A2 - C2 * A1) / (A1 * B2 - A2 * B1)};
}

bool sweepline::event_t::operator < (const event_t &e) const {
    if(std::fabs(p.x - e.p.x) > geometry::EPS)
        return p.x < e.p.x;
    else if(std::fabs(p.y - e.p.y) > geometry::EPS)
        return p.y < e.p.y;
    else
        return seg_id < e.seg_id;
}

std::vector<sweepline::intersection_t> sweepline::find_intersections(
        const std::vector<geometry::segment_t> &line_segments,
        bool verbose,
        bool enable_color
) {

    return sweepline::solver(line_segments, verbose, enable_color).solve();
}

geometry::float_t sweepline::sweeplineX;

sweepline::solver::solver(const std::vector<geometry::segment_t> &line_segments, bool verbose, bool enable_color)
        : line_segments(line_segments), verbose(verbose) {
}

std::vector<sweepline::intersection_t> sweepline::solver::solve() {
    // initialize the sweepline to -inf
    sweepline::sweeplineX = -std::numeric_limits<geometry::float_t>::max();

    // initialize the event_queue by inserting the end points of the line segments
    // and populate vertical_segs with vertical segments
    init_event_queue();

    // sort vertical segments by x then y
    std::sort(vertical_segs.begin(), vertical_segs.end(),
              [](const geometry::segment_t &a, const geometry::segment_t &b) {
                  return a.p.x == b.p.x? a.p.y < b.p.y : a.p.x < b.p.x;
              }
    );

    // find intersections between pairs of vertical line segments
    find_vertical_vertical_intersections();

    while(!event_queue.empty()) {
        sweepline::event_t top = *event_queue.begin();
        event_queue.erase(event_queue.begin());

        if(top.p.x < sweepline::sweeplineX) {
            continue;
        }

        // find intersections between a vertical line segment and one or more non-vertical segments
        find_vertical_nonvertical_intersections(top.p.x);

        // move sweepline to x coordinate of event being processed
        sweepline::sweeplineX = top.p.x;

        // get the active segments with an event at the point currently being processed
        // returns three arrays of active segment indices corresponding to event_t::type
        auto active_segs = get_active_segs(top);

        // remove all end points, insert all begin points and reorder the interior points
        update_segment_ordering(active_segs);

        // if no segments were newly inserted, the immediate left and right neighbours
        // of the deleted set of segments become adjacent candidates for intersection
        if(active_segs[sweepline::event_t::type::begin].size()
           + active_segs[sweepline::event_t::type::interior].size() == 0)
            handle_no_newly_inserted(top.p);
        else
            handle_extremes_of_newly_inserted();
        // else the left and right extremes among the set of newly inserted segments
        // must be checked for intersection with the immediate left and right neighbours respectively

        // reset the sweepline to the x coordinate of event being processed
        sweepline::sweeplineX = top.p.x;

        // finally report the union of all active segments as an intersection if there are two or more of them
        if(active_segs[0].size() + active_segs[1].size() + active_segs[2].size() > 1)
            report_intersection(top.p, std::move(active_segs));
    }

    merge_intersection_points();

    if(verbose)
        std::cerr << std::endl;

    return result;
}

void sweepline::solver::init_event_queue() {
    for(size_t i = 0; i < line_segments.size(); i++) {
        const auto &p = line_segments[i].p;
        const auto &q = line_segments[i].q;

        if(std::fabs(p.x - q.x) < geometry::EPS) {
            // handle (vertical) segments with same slope as sweepline separately
            vertical_segs.emplace_back(line_segments[i]);
        } else {
            // insert the begin and end points of each segment in the event queue
            event_queue.insert({ p, sweepline::event_t::type::begin, i});
            event_queue.insert({ q, sweepline::event_t::type::end,   i});
        }
    }
}

void sweepline::solver::find_vertical_vertical_intersections() {
    for(size_t i = 0; i + 1 < vertical_segs.size(); i++) {
        if(vertical_segs[i].q == vertical_segs[i + 1].p) {
            sweepline::intersection_t it {
                    vertical_segs[i].q,
                    std::vector<size_t> {
                            vertical_segs[i].seg_id,
                            vertical_segs[i + 1].seg_id
                    }
            };

            result.emplace_back(it);
        }
    }
}

void sweepline::solver::find_vertical_nonvertical_intersections(geometry::float_t max_vsegx) {
    while(vert_idx < vertical_segs.size()
          and vertical_segs[vert_idx].p.x < sweepline::sweeplineX - geometry::EPS)
        vert_idx++;

    while(vert_idx < vertical_segs.size()
          and vertical_segs[vert_idx].p.x <= max_vsegx + geometry::EPS) {

        auto &vseg = vertical_segs[vert_idx];
        sweepline::sweeplineX = vseg.p.x;

        auto itr = seg_ordering.lower_bound(geometry::segment_t{ vseg.p, vseg.p, 0 });

        while(itr != seg_ordering.end()) {
            geometry::float_t it_y = itr->eval_y(sweepline::sweeplineX);

            if(it_y > vseg.q.y + geometry::EPS)
                break;

            sweepline::intersection_t it {
                    geometry::point_t{ sweepline::sweeplineX, it_y },
                    std::vector<size_t>{ itr->seg_id, vseg.seg_id }
            };

            result.emplace_back(it);

            ++itr;
        }

        vert_idx++;
    }
}

std::array<std::vector<size_t>, 3> sweepline::solver::get_active_segs(sweepline::event_t top) {
    // array of all segments with an event at the point currently being processed
    //   active[event_t::type::begin]    -> list of segments which begin at this point
    //   active[event_t::type::interior] -> list of segments which intersect with some other segment at this point
    //   active[event_t::type::end]      -> list of segments which end at this point
    std::array<std::vector<size_t>, 3> active;
    active[top.tp].push_back(top.seg_id);

    // get all segments with an event at top.p and add them to one of the above
    while(!event_queue.empty() and event_queue.begin()->p == top.p) {
        sweepline::event_t nxt_top = *event_queue.begin();
        event_queue.erase(event_queue.begin());
        active[nxt_top.tp].push_back(nxt_top.seg_id);
    }

    return active;
}

void sweepline::solver::update_segment_ordering(const std::array<std::vector<size_t>, 3> &active) {
    // remove all end event segments
    for(size_t idx: active[sweepline::event_t::type::end])
        seg_ordering.erase(line_segments[idx]);

    // remove all interior event segments
    for(size_t idx: active[sweepline::event_t::type::interior])
        seg_ordering.erase(line_segments[idx]);

    // increment the sweepline by a very small amount, just past the intersection point
    sweepline::sweeplineX += geometry::EPS_INC;

    max_y = -std::numeric_limits<geometry::float_t>::max();
    min_y = std::numeric_limits<geometry::float_t>::max();

    // insert all begin type events
    for(int idx: active[sweepline::event_t::type::begin]) {
        min_y = std::min(min_y, line_segments[idx].eval_y(sweepline::sweeplineX));
        max_y = std::max(max_y, line_segments[idx].eval_y(sweepline::sweeplineX));
        seg_ordering.insert(line_segments[idx]);
    }

    // re-insert all interior type events (so that ordering is updated)
    for(int idx: active[sweepline::event_t::type::interior]) {
        min_y = std::min(min_y, line_segments[idx].eval_y(sweepline::sweeplineX));
        max_y = std::max(max_y, line_segments[idx].eval_y(sweepline::sweeplineX));
        seg_ordering.insert(line_segments[idx]);
    }
}

void sweepline::solver::handle_no_newly_inserted(geometry::point_t cur) {
    auto b_right = seg_ordering.lower_bound(geometry::segment_t{ cur, cur, 0 });
    if(b_right != seg_ordering.end() and b_right != seg_ordering.begin()) {
        auto b_left = b_right;
        --b_left;
        if(geometry::is_intersecting(*b_left, *b_right)) {
            geometry::point_t pt = geometry::intersection_point(*b_left, *b_right);
            sweepline::event_t::type tp1 = b_left->p == pt? sweepline::event_t::type::begin : sweepline::event_t::type::interior;
            sweepline::event_t::type tp2 = b_right->p == pt? sweepline::event_t::type::begin : sweepline::event_t::type::interior;

            event_queue.insert(sweepline::event_t{ pt, tp1, b_left->seg_id  });
            event_queue.insert(sweepline::event_t{ pt, tp2, b_right->seg_id });
        }
    }
}

void sweepline::solver::handle_extremes_of_newly_inserted() {
    geometry::point_t left { sweepline::sweeplineX, min_y - 2 * geometry::EPS }, right { sweepline::sweeplineX, max_y + 2 * geometry::EPS };
    auto b_right = seg_ordering.lower_bound(geometry::segment_t{ right, right, 0 });
    auto s_left  = seg_ordering.lower_bound(geometry::segment_t{ left,  left,  0 });

    // check for candidate intersection at the right extreme
    if(b_right != seg_ordering.end() and b_right != seg_ordering.begin()) {
        auto s_right = b_right;
        --s_right;

        if(is_intersecting(*s_right, *b_right)) {
            geometry::point_t pt = intersection_point(*s_right, *b_right);
            sweepline::event_t::type tp1 = s_right->p == pt? sweepline::event_t::type::begin : sweepline::event_t::type::interior;
            sweepline::event_t::type tp2 = b_right->p == pt? sweepline::event_t::type::begin : sweepline::event_t::type::interior;

            event_queue.insert(sweepline::event_t{ pt, tp1, s_right->seg_id });
            event_queue.insert(sweepline::event_t{ pt, tp2, b_right->seg_id });
        }
    }

    // check for candidate intersection at the left extreme
    if(s_left != seg_ordering.end() and s_left != seg_ordering.begin()) {
        auto b_left = s_left;
        --b_left;

        if(is_intersecting(*b_left, *s_left)) {
            geometry::point_t pt = intersection_point(*b_left, *s_left);
            sweepline::event_t::type tp1 = b_left->p == pt? sweepline::event_t::type::begin : sweepline::event_t::type::interior;
            sweepline::event_t::type tp2 = s_left->p == pt? sweepline::event_t::type::begin : sweepline::event_t::type::interior;

            event_queue.insert(sweepline::event_t{ pt, tp1, b_left->seg_id });
            event_queue.insert(sweepline::event_t{ pt, tp2, s_left->seg_id });
        }
    }
}

void sweepline::solver::report_intersection(geometry::point_t cur, std::array<std::vector<size_t>, 3> &&active) {
    active[0].insert(active[0].end(), active[1].begin(), active[1].end());
    active[0].insert(active[0].end(), active[2].begin(), active[2].end());

    result.emplace_back(sweepline::intersection_t{ cur, std::move(active[0]) });
}

void sweepline::solver::merge_intersection_points() {
    std::sort(result.begin(), result.end(),
              [](const sweepline::intersection_t &a, const sweepline::intersection_t &b) {
                  return a.pt.x == b.pt.x? a.pt.y < b.pt.y : a.pt.x < b.pt.x;
              }
    );

    std::vector<sweepline::intersection_t> merged;
    for(size_t i = 0, j; i < result.size(); i = j) {
        std::vector<size_t> indices;
        for(j = i; j < result.size() and result[j].pt == result[i].pt; j++)
            indices.insert(indices.end(), result[j].segments.begin(), result[j].segments.end());

        std::sort(indices.begin(), indices.end());
        indices.erase(std::unique(indices.begin(), indices.end()), indices.end());

        merged.emplace_back(
                sweepline::intersection_t {
                        result[i].pt, indices
                }
        );
    }

    result = std::move(merged);
}

int main() {
    freopen("A.in", "r", stdin);
    size_t n;       // number of input segments
    std::cin >> n;

    std::vector<geometry::segment_t> segments;
    segments.reserve(n);

    for(size_t i = 0; i < n; i++) {
        float_t x1, y1, x2, y2;
        std::cin >> x1 >> y1 >> x2 >> y2;

        if(std::make_pair(x1, y1) > std::make_pair(x2, y2))
            std::swap(x1, x2), std::swap(y1, y2);

        segments.emplace_back(geometry::segment_t{
                geometry::point_t{ x1, y1 },
                geometry::point_t{ x2, y2 },
                i });
    }

    std::vector<sweepline::intersection_t> result = sweepline::find_intersections(segments, false, false);
    cout << result.size() << endl;
}
