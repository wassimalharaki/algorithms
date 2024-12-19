// O(n^2log(n))
//import math
//
//# Utility function to compute the angle between two points
//def angle_between(p1, p2):
//# Compute the angle in radians between point p1 and p2
//return math.atan2(p2[1] - p1[1], p2[0] - p1[0])
//
//# Function to solve the problem using Radial Sweep
//def min_deletions_to_separate_points(red_points, blue_points):
//N = len(red_points)
//M = len(blue_points)
//
//# Try using each red point as the reference point
//min_deletions = float('inf')
//
//for ref_red in red_points:
//events = []
//
//# Calculate angles for all points relative to the current reference red point
//for red in red_points:
//if red != ref_red:
//events.append((angle_between(ref_red, red), 'red'))
//
//for blue in blue_points:
//events.append((angle_between(ref_red, blue), 'blue'))
//
//# Sort events by angle
//events.sort()
//
//# Initialize counts of points on each side of the line
//red_on_right = N - 1  # All reds except the reference are initially considered on the right
//        blue_on_right = M  # All blues are initially considered on the right
//
//# Sweep through all events (angles), and update counts
//for angle, color in events:
//if color == 'red':
//red_on_right -= 1  # One red point moves to the left
//else:
//blue_on_right -= 1  # One blue point moves to the left
//
//# Calculate how many deletions we need to achieve separation
//        deletions = red_on_right + (M - blue_on_right)  # Red deletions + Blue deletions
//min_deletions = min(min_deletions, deletions)
//
//return min_deletions
