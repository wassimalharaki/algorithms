// O(n sqrt(n))
#define double long double
#define px second
#define py first
typedef pair<long long, long long> pairll;
const int MAX = 1e5;
pairll pnts [MAX];
int compare(pairll a, pairll b)
{
    return a.px<b.px;
}
double closest_pair(pairll pnts[],int n) {
    sort(pnts,pnts+n,compare);
    double best=INF;
    set<pairll> box;
    box.insert(pnts[0]);
    int left = 0;
    for (int i=1;i<n;++i)
    {
        while (left<i && pnts[i].px-pnts[left].px > best)
            box.erase(pnts[left++]);
        for(typeof(box.begin()) it=box.lower_bound(make_pair(pnts[i].py-best, pnts[i].px-best));it!=box.end() && pnts[i].py+best>=it->py;it++)
            best = min((double)best, (double)sqrt(pow(pnts[i].py - it->py, 2.0)+pow(pnts[i].px - it->px, 2.0)));
        box.insert(pnts[i]);
    }
    return best;
}
