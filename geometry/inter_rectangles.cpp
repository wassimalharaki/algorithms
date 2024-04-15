#include <bits/stdc++.h>
using namespace std;

//left, down, right, up
//if either val is < 0, no intersection
//if both vals are = 0, single point intersection
//if either is 0, side intersection
//else inter_area > 0
// O(1)
int inter_area(vector<int>& r1, vector<int>& r2) {
	return max(0, min(r1[2], r2[2]) - max(r1[0], r2[0])) *
	       max(0, min(r1[3], r2[3]) - max(r1[1], r2[1]));
}