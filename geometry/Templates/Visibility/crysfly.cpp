// source https://qoj.ac/submission/626090
// what is matter? never mind.
#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")
//#pragma GCC target("sse,sse2,sse3,sse4,popcnt,abm,mmx,avx,avx2")
#include<bits/stdc++.h>
#define For(i,a,b) for(int i=(a);i<=(b);++i)
#define Rep(i,a,b) for(int i=(a);i>=(b);--i)
#define ll long long
#define ull unsigned long long
//#define int long long
#define SZ(x) ((int)((x).size()))
#define ALL(x) (x).begin(),(x).end()
using namespace std;
inline int read()
{
    char c=getchar();int x=0;bool f=0;
    for(;!isdigit(c);c=getchar())f^=!(c^45);
    for(;isdigit(c);c=getchar())x=(x<<1)+(x<<3)+(c^48);
    return f?-x:x;
}

#define fi first
#define se second
#define pb push_back
#define mkp make_pair
typedef pair<int,int>pii;
typedef vector<int>vi;

typedef long double db;
const db eps=1e-8,pi=3.14159265358979323846;
int sgn(db x){return x<-eps?-1:x>eps;}
int cmp(db a,db b){return sgn(a-b);}

struct P{
    db x,y;
    P(db x=0,db y=0):x(x),y(y){}
    P&operator +=(P o){return x+=o.x,y+=o.y,*this;}
    P&operator -=(P o){return x-=o.x,y-=o.y,*this;}
    P&operator *=(db o){return x*=o,y*=o,*this;}
    P&operator /=(db o){return x/=o,y/=o,*this;}
    friend P operator +(P a,P b){return a+=b;}
    friend P operator -(P a,P b){return a-=b;}
    friend P operator *(P a,db b){return a*=b;}
    friend P operator /(P a,db b){return a/=b;}
    friend bool operator <(P a,P b){return fabs(a.x-b.x)<eps?a.y<b.y:a.x<b.x;}
    friend bool operator ==(P a,P b){return cmp(a.x,b.x)==0 && cmp(a.y,b.y)==0;}
    friend bool operator !=(P a,P b){return !(a==b);}
    friend db operator %(P a,P b){return a.x*b.x+a.y*b.y;} // dot
    friend db operator *(P a,P b){return a.x*b.y-a.y*b.x;} // cross

    P rot(db o){
        db s=sin(o),c=cos(o),xx=x*c-y*s,yy=x*s+y*c;
        x=xx,y=yy;return *this;
    }
    P rot90(){swap(x,y),x=-x;return *this;}
    db ang(){return atan2(y,x);}
    db len(){return sqrt(x*x+y*y);}
    db len2(){return x*x+y*y;}

    int half(){return sgn(y)==1||(sgn(y)==0&&sgn(x)>=0);}
    P unit(){return ((*this))/len();}

    void read(){cin>>x>>y;}
    void out(){cout<<"("<<x<<","<<y<<")"<<endl;}
};
bool cmp_dir(P a,P b){
    if(a.half()!=b.half())return a.half()<b.half();
    return sgn(a*b)>0;
}

db dis(P a,P b){return (a-b).len();}
db cross(P a,P b,P c){
    // (a->b)*(a->c)
    return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
}
int cmp3(P a,P b,P c){
    return sgn(cross(a,b,c));
}

bool in_tri(P a,P b,P c,P p){
    if(cmp3(a,b,c)<0) swap(b,c);
    return cmp3(a,b,p)>=0 && cmp3(b,c,p)>=0 && cmp3(c,a,p)>=0;
}
db area_tri(P a,P b,P c){
    return fabs(cross(a,b,c))/2;
}

bool paral(P p1,P p2,P q1,P q2){
    // is parallel
    return sgn((p2-p1)*(q2-q1))==0;
}
P interl(P p1,P p2,P q1,P q2){
    // intersect point
    db s1=cross(q1,q2,p1),s2=-cross(q1,q2,p2);
    return (p1*s2+p2*s1)/(s1+s2);
}

bool inter(db l1,db r1,db l2,db r2){
    if(l1>r1)swap(l1,r1); if(l2>r2)swap(l2,r2);
    return !(cmp(r1,l2)==-1||cmp(r2,l1)==-1);
}
bool ismid(db a,db m,db b){
    return sgn(a-m)==0||sgn(b-m)==0||((a<m)!=(b<m));
}
bool ismid(P a,P m,P b){
    return ismid(a.x,m.x,b.x)&&ismid(a.y,m.y,b.y);
}

bool isseg(P p1,P p2,P q1,P q2){
    return inter(p1.x,p2.x,q1.x,q2.x) && inter(p1.y,p2.y,q1.y,q2.y) &&
           cmp3(p1,p2,q1)*cmp3(p1,p2,q2)<=0 && cmp3(q1,q2,p1)*cmp3(q1,q2,p2)<=0;
}
bool isseg_strict(P p1,P p2,P q1,P q2){
    return cmp3(p1,p2,q1)*cmp3(p1,p2,q2)<0 && cmp3(q1,q2,p1)*cmp3(q1,q2,p2)<0;
}

struct L{
    P a,b;
    L(P aa,P bb){a=aa,b=bb;}
    bool in(P p){return sgn((b-a)*(p-a))>0;}
    int in_sgn(P p){return sgn((b-a)*(p-a));}
    P dir(){return b-a;}
    bool onl(P p){
        return cmp3(a,b,p)==0;
    }
    bool onseg(P p){
        return onl(p)&&ismid(a,p,b);
    }
    bool onseg_strict(P p){
        return onl(p)&&sgn((p-a)%(a-b))*sgn((p-b)%(a-b))<0;
    }
    void out(){cout<<"("<<a.x<<","<<a.y<<")---("<<b.x<<","<<b.y<<")\n";}
};

bool isseg(L a,L b){
    return isseg(a.a,a.b,b.a,b.b);
}
bool paral(L a,L b){
    // is parallel
    return paral(a.a,a.b,b.a,b.b);
}
P operator &(L a,L b){
    return interl(a.a,a.b,b.a,b.b);
}
bool samedir(L a,L b){
    return paral(a,b) && sgn(a.dir()%b.dir())==1;
}
bool operator <(L a,L b){
    if(samedir(a,b)) return b.in(a.a);
    return cmp_dir(a.dir(),b.dir());
}

P proj(L a,P b){
    P d=a.dir();
    return a.a+d*((d%(b-a.a))/d.len2());
}
P reflect(L a,P b){
    return proj(a,b)*2-b;
}
db dis(L a,P b){
    db s=abs((b-a.a)*(b-a.b));
    return s/dis(a.a,a.b);
}
db dis_seg(L a,P b){
    if(a.a==a.b) return	dis(a.a,b);
    P h=proj(a,b);
    if(ismid(a.a,h,a.b)) return dis(h,b);
    return min(dis(a.a,b),dis(a.b,b));
}

db rad(P a,P b){
    return atan2l(a*b,a%b);
}

// polygon
db area(vector<P>a){
    db res=0;
    For(i,0,(int)a.size()-1)res+=a[i]*a[(i+1)%a.size()];
    return res/2;
}
int contain(vector<P>a,P p){
    int n=a.size(),res=0;
    For(i,0,n-1){
        P u=a[i],v=a[(i+1)%n];
        if(L(u,v).onseg(p))return 1;
        if(cmp(u.y,v.y)<=0)swap(u,v);
        if(cmp(p.y,u.y)>0 || cmp(p.y,v.y)<=0)continue;
        res^=cmp3(p,u,v)>0;
    }
    return res*2;
}
vector<P>cut(vector<P>&a,P q1,P q2){
    vector<P>b; int n=a.size();
    For(i,0,n-1){
        P p1=a[i],p2=a[(i+1)%n];
        int d1=cmp3(q1,q2,p1),d2=cmp3(q1,q2,p2);
        if(d1>=0)b.pb(p1);
        if(d1*d2<0)b.pb(interl(p1,p2,q1,q2));
    }
    return b;
}
vector<P>cut(vector<P>&a,L l){
    return cut(a,l.a,l.b);
}

vector<P>convex(vector<P>a){
    int n=a.size(),m=0; if(n<=1)return a;
    sort(a.begin(),a.end());
    vector<P>st(n*2); int tp=0;
    For(i,0,n-1){
        while(tp>1 && cmp3(st[tp-2],st[tp-1],a[i])<=0)--tp;
        st[tp++]=a[i];
    }
    int t=tp;
    Rep(i,n-2,0){
        while(tp>t && cmp3(st[tp-2],st[tp-1],a[i])<=0)--tp;
        st[tp++]=a[i];
    }
    st.resize(tp-1);
    return st;
}
db diam(vector<P>a){
    int n=a.size();
    if(n<=1)return 0;
    int ii=0,jj=0;
    For(k,1,n-1){
        if(a[k]<a[ii])ii=k;
        if(a[jj]<a[k])jj=k;
    }
    int i=ii,j=jj;
    db res=dis(a[i],a[j]);
    do{
        if((a[(i+1)%n]-a[i])*(a[(j+1)%n]-a[j])>=0) (++j)%=n;
        else (++i)%=n;
        res=max(res,dis(a[i],a[j]));
    }while(i!=ii||j!=jj);
    return res;
}

bool check(L a,L b,L c){
    return c.in(a&b);
}
int checksgn(L a,L b,L c){
    return c.in_sgn(a&b);
}
vector<P>hpis(vector<L>&l){
    sort(l.begin(),l.end());
    deque<L>q;
    For(i,0,(int)l.size()-1){
        if(i&&samedir(l[i],l[i-1]))continue;
        while(q.size()>1 && !check(q[q.size()-2],q[q.size()-1],l[i]))q.pop_back();
        while(q.size()>1 && !check(q[1],q[0],l[i]))q.pop_front();
        q.pb(l[i]);
    }
    while(q.size()>2 && !check(q[q.size()-2],q[q.size()-1],q[0]))q.pop_back();
    while(q.size()>2 && !check(q[1],q[0],q[q.size()-1]))q.pop_front();
    vector<P>res;
    For(i,0,(int)q.size()-1) res.pb(q[i]&q[(i+1)%q.size()]);
    return res;
}

mt19937_64 rnd(64);

vector<L>cut(vector<L>&a,L l){
    vector<L>b; int n=a.size();
    For(i,0,n-1){
        L a1=a[i],a2=a[(i+1)%n],a3=a[(i+2)%n];
        int d1=checksgn(a1,a2,l),d2=checksgn(a2,a3,l);
        if(d1>0 || d2>0 || (d1==0&&d2==0)) b.pb(a2);
        if(d1>=0 && d2<0) b.pb(l);
    }
    return b;
}

#define maxn 500005
#define inf 0x3f3f3f3f

int n;
P a[255];
db f[105][105],g[105];

bool inangle(P a,P b,P c){
    if(sgn(a*b)>0) return sgn(a*c)>=0 && sgn(c*b)>=0;
    return !(sgn(b*c)>0 && sgn(c*a)>0);
}
bool chk(P s,P t,bool half=0){
    if(s==t)return 1;
    vi cnt(2);
    For(i,0,n-1){
        P a=::a[i],b=::a[(i+1)%n];
        if(isseg_strict(s,t,a,b)) return 0;
        if(L(s,t).onseg(a)){
            P c=::a[(i+n-1)%n];
            if(s!=a && !inangle(s-a,b-a,c-a)) return 0;
            if(t!=a && !inangle(t-a,b-a,c-a)) return 0;
            if(half && s!=a && t!=a){
                int d=sgn((c-a)*(t-s));
                if(d) cnt[d>0]++;
                d=sgn((b-a)*(t-s));
                if(d) cnt[d>0]++;
                if(cnt[1] && cnt[0]) return 0;
            }
        }else{
            if(s!=a && s!=b && L(a,b).onseg(s) && sgn((t-s)*(b-a))>0) return 0;
            if(t!=a && t!=b && L(a,b).onseg(t) && sgn((s-t)*(b-a))>0) return 0;
        }
    }
    return 1;
}

bool vis[105];

signed main()
{
    freopen("A.in", "r", stdin);
    n=read();
    For(i,0,n+1)a[i].x=read(),a[i].y=read();
    if(chk(a[n],a[n+1],1)) {
        puts("0.0000");
        exit(0);
    }
    For(i,0,n+1){
        f[i][i]=0;
        For(j,i+1,n+1) {
            if(chk(a[i],a[j])) f[i][j]=dis(a[i],a[j]);
            else f[i][j]=4e18;
            f[j][i]=f[i][j];
        }
    }
    For(k,0,n+1) For(i,0,n+1) For(j,0,n+1) f[i][j]=min(f[i][j],f[i][k]+f[k][j]);
    For(i,0,n+1) g[i]=4e18;
    For(i,0,n-1){
        P d=a[i]-a[n+1];
        P dd=d; dd.rot90();
        L t1=L(a[n+1],a[i]);
        For(j,0,n+1){
            L t2=L(a[j],a[j]+dd);
            if(!paral(t1,t2)){
                P p=t1&t2;
                if(chk(a[j],p,0) && chk(p,a[n+1],1)) {
                    g[j]=min(g[j],dis(p,a[j]));
                }
            }
        }
        For(j,0,n-1){
            L t2=L(a[j],a[(j+1)%n]);
            if(!paral(t1,t2)){
                P p=t1&t2;
                if(chk(p,a[n+1],1)) {
                    //		cout<<"go "<<j<<" "<<(j+1)%n<<"\n";
                    For(k,0,n+1)
                        if(chk(a[k],p,0)){
                            g[k]=min(g[k],dis(p,a[k]));
                        }
                }
            }
        }
    }
    For(i,0,n+1) cout<<f[i][n]<<" "<<g[i]<<"\n";
    db res=4e18;
    For(i,0,n+1) res=min(res,g[i]+f[i][n]);
    printf("%.12lf\n",(double)res);
    return 0;
}