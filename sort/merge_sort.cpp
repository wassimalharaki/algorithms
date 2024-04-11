#include <bits/stdc++.h>
using namespace std;

void merge(vector<int>& nums, int left, int mid, int right) {
    int ls = mid - left + 1, rs = right - mid;
    vector<int> subLeft(ls), subRight(rs);

    vector<int>::iterator begin = nums.begin();
    copy(begin + left, begin + left + ls, subLeft.begin());
    copy(begin + mid + 1, begin + mid + 1 + rs, subRight.begin());
    
    int i = 0, j = 0, k = left;
    while (i < ls && j < rs)
        subLeft[i] < subRight[j] ?
            nums[k++] = subLeft[i++]:
            nums[k++] = subRight[j++];
    
    while (i < ls)
        nums[k++] = subLeft[i++];
    while (j < rs)
        nums[k++] = subRight[j++];
}

void mergeSort(vector<int>& nums, int left, int right) { 
    if (left >= right)
        return;

    int mid = (left + right) / 2;
    mergeSort(nums, left, mid); 
    mergeSort(nums, mid + 1, right); 
    merge(nums, left, mid, right);
}