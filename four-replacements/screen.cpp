// Exact necessary finite hitting-set test on single-deletion time windows.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>
using Z=std::int64_t;
using Wide=__int128;
using Mask=std::uint64_t;
using Row=std::vector<int>;
struct Point{Z p,q;}; // scaled time z=n*t
struct Interval{Point a,b;};
Point frac(Z p,Z q){Z g=std::gcd(p,q); return {p/g,q/g};}
Point blend(Point a,Point b,int t,int d){
    return frac((d-t)*a.p*b.q+t*b.p*a.q,d*a.q*b.q);
}
int inverse(int a,int m){
    int r=m,s=a,u=0,v=1;
    while(s){int q=r/s;int t=r-q*s;r=s;s=t;t=u-q*v;u=v;v=t;}
    if(r!=1)throw std::logic_error("nonunit");
    return (u%m+m)%m;
}
struct Family{
    int r=0;
    std::vector<Row> rows;
    std::vector<Mask> signature;
    Mask all=0;
};
// Equal-length windows at several primitive centres. Their covered lengths can
// be added even when the windows overlap, since every window must be covered.
bool bundle(int n,int r,int low,int high,const Row& centres,Row& answer){
    constexpr int d=24;
    if(centres.empty())return false;
    Z cap=8*Z(n)*r*d/((n-8)*(high-low));
    if(cap>80*n || cap>80000)return false;
    Z period=Z(n)*r*d, radius=r*d;
    for(int v=n;v<=cap;++v){
        Z width=Z(v)*(high-low), total=0;
        for(int p:centres){
            Z start=Z(n)*d*((Z(v)*p)%r)+Z(v)*low, end=start+width;
            for(Z k=(start+radius)/period;(k*period-radius)<end;++k){
                Z a=std::max(start,k*period-radius), b=std::min(end,k*period+radius);
                if(b>a)total+=b-a;
            }
        }
        if(4*total>=Z(centres.size())*width)answer.push_back(v);
    }
    return true;
}
// All speeds which cover at least one fourth of this interval's length.
// In coordinates w=v*q_a*q_b*z the interval and every band have integer endpoints.
Row quarter(int n,Interval c){
    Point a=c.a,b=c.b;
    Z den=a.q*b.q,len=b.p*a.q-a.p*b.q;
    if(len<=0)throw std::logic_error("empty probe");
    Z cap=8*Z(n)*den/((n-8)*len);
    if(cap>80*n || cap>80000)return {}; // skip entire row (a relaxation)
    Row row;
    for(int v=n;v<=cap;++v){
        Z left=Z(v)*a.p*b.q, right=Z(v)*b.p*a.q;
        Z j=(left+den)/(n*den);
        Z mass=0;
        for(; (n*j-1)*den<right; ++j){
            Z lo=std::max(left,(n*j-1)*den);
            Z hi=std::min(right,(n*j+1)*den);
            if(hi>lo)mass+=hi-lo;
        }
        if(4*mass>=right-left)row.push_back(v);
    }
    return row;
}
Family family(int n,int r,bool budget){
    Family f; f.r=r;
    auto add=[&](Interval c){
        if(f.rows.size()>=60)return;
        Row row=quarter(n,c);
        if(!row.empty() && std::find(f.rows.begin(),f.rows.end(),row)==f.rows.end())
            f.rows.push_back(std::move(row));
    };
    if(r==1){
        add({{1,2},{1,1}});
        add({{1,2},{3,4}});
        add({{3,4},{1,1}});
    }else{
        int s=n-r,used=0;
        constexpr int slices[][2]={{0,12},{6,12},{12,18},{0,18},{4,12},{8,16},
            {12,16},{12,20},{6,18},{16,20},{0,8},{3,9},{8,20},{14,22},{12,24}};
        for(auto& slice:slices){
            int low=slice[0],high=slice[1];
            if(2*r<n && low<12)continue;
            Row centres;
            for(int F=s;F<n && high*F<=24*s;++F)if(std::gcd(F,r)==1){
                centres.push_back(r-inverse(F%r,r));
                if(centres.size()>=24)break;
            }
            if(centres.size()<2)continue;
            Row row;
            if(bundle(n,r,low,high,centres,row)){
                if(row.empty()){f.rows={Row{}}; f.all=1;return f;}
                f.rows.push_back(std::move(row));
            }
        }
        for(int F=s;F<n && used<24;++F){
            if(std::gcd(r,F)!=1)continue;
            int p=r-inverse(F%r,r);
            Point start=2*r>=n ? frac(Z(n)*p,r) : frac(2*Z(n)*p+1,2*r);
            Point end=frac(Z(n)*((Z(F)*p+1)/r)-1,F);
            if(end.p*start.q<=start.p*end.q)continue;
            ++used;
            add({start,end});
            // Each positive subinterval also has to be covered.
            if(used<=12){
                Point mid=blend(start,end,1,2);
                add({start,mid});add({mid,end});
            }
        }
    }
    if(2*r>=n){
        Row multiples;
        int maxv=std::min(80*n,80000);
        if(!budget || Z(n)*r/(2*(n-r))>maxv)multiples.push_back(0); // unrestricted tail
        for(int m=2;m*r<=maxv && (!budget || 2*m*(n-r)<=n);++m)multiples.push_back(m*r);
        f.rows.push_back(std::move(multiples));
    }
    std::sort(f.rows.begin(),f.rows.end(),[](const Row&a,const Row&b){return a.size()<b.size();});
    int maxv=0;
    for(const Row& row:f.rows)if(!row.empty())maxv=std::max(maxv,row.back());
    f.signature.assign(maxv+1,0);
    for(int i=0;i<int(f.rows.size());++i)
        for(int v:f.rows[i])f.signature[v]|=Mask(1)<<i;
    f.all=(Mask(1)<<f.rows.size())-1;
    return f;
}
using Prefix=std::array<const Family*,4>;
using Remaining=std::array<Mask,4>;
struct Solver{
    Prefix families;
    int count,marked,n;
    bool allowed(int v)const{
        if(!v || !marked || 2*marked<n || v%marked)return true;
        return 2*Z(v/marked)*(n-marked)>n;
    }
    bool search(Remaining need,int left)const{
        const Row* shortest=nullptr;
        for(int i=0;i<count;++i){
            if(!need[i])continue;
            int j=__builtin_ctzll(need[i]);
            const Row& row=families[i]->rows[j];
            if(!shortest || row.size()<shortest->size())shortest=&row;
        }
        if(!shortest)return true;
        if(!left)return false;
        for(int v:*shortest)if(allowed(v)){
            Remaining next=need;
            for(int i=0;i<count;++i){
                const auto& m=families[i]->signature;
                if(v<int(m.size()))next[i]&=~m[v];
            }
            if(search(next,left-1))return true;
        }
        return false;
    }
    bool feasible()const{
        Remaining all{};
        for(int i=0;i<count;++i)all[i]=families[i]->all;
        return search(all,4);
    }
};
int main(int argc,char**argv){
    if(argc!=4)return 2;
    int lo=std::stoi(argv[2]),hi=std::stoi(argv[3]);
    if(lo<9||hi>2500||hi<lo)return 2;
    std::ifstream in(argv[1]); if(!in)return 2;
    int n,good,bad;
    while(in>>n>>good>>bad){
        Row rs(good+bad);
        for(int& r:rs)if(!(in>>r))return 2;
        if(n<lo||n>hi)continue;
        std::vector<Family> data(rs.size());
        auto get=[&](int i){
            if(!data[i].r)data[i]=family(n,rs[i],i<good);
            return &data[i];
        };
        Z tried=0,survived=0;
        // Root at one deletion which lacks a small multiple, or lies in the lower half.
        for(int q=good;q<int(rs.size());++q){
            Solver solve{{get(q),nullptr,nullptr,nullptr},1,rs[q],n};
            if(!solve.feasible())continue;
            for(int a=0;a<int(rs.size());++a)if(a!=q){
                solve.families[1]=get(a); solve.count=2;
                if(!solve.feasible())continue;
                for(int b=a+1;b<int(rs.size());++b)if(b!=q){
                    solve.families[2]=get(b);solve.count=3;
                    if(!solve.feasible())continue;
                    for(int c=b+1;c<int(rs.size());++c)if(c!=q){
                        solve.families[3]=get(c); solve.count=4;++tried;
                        if(solve.feasible()){
                            ++survived;
                            std::cout<<"Q "<<n<<' '<<rs[q]<<' '<<rs[a]<<' '<<rs[b]<<' '<<rs[c]<<'\n';
                        }
                    }
                }
            }
        }
        std::cout<<"N "<<n<<' '<<tried<<' '<<survived<<std::endl;
    }
    if(!in.eof())return 2;
}
