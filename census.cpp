// cl /std:c++17 /O2 /EHsc census.cpp
// g++ -std=c++17 -O3 census.cpp -o census
// Enumerate every two-entry modification with LR <= 1/n for 8 <= n <= 300.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using Int = std::int64_t;
struct Q { Int p, d; }; // x = n*t, 0 <= x <= n/2
bool operator<(Q a, Q b) { return a.p*b.d < b.p*a.d; }
bool equal(Q a, Q b) { return !(a<b) && !(b<a); }
struct Interval { Q a, b; };
Q length(Interval i) { return {i.b.p*i.a.d-i.a.p*i.b.d, i.a.d*i.b.d}; }
struct Event { Q x; int v, change; };
using List = std::vector<Interval>;

void append(List& out, Interval i) {
    if (!(i.a<i.b)) return;
    if (!out.empty() && equal(out.back().b,i.a)) out.back().b=i.b;
    else out.push_back(i);
}

// Between consecutive walls record intervals with one or two bad original
// speeds. Dirichlet's pigeonhole argument precludes a cell with zero.
std::map<std::pair<int,int>,List> cells(int n) {
    std::vector<Event> events;
    std::set<int> active;
    Q stop{n,2};
    for(int v=1;v<n;++v){
        active.insert(v); // zero is in B_v
        for(int j=0;j<=v/2+1;++j){
            Q a{Int(n)*j-1,v}, b{Int(n)*j+1,v};
            if(a.p>0 && a<stop) events.push_back({a,v,1});
            if(b<stop) events.push_back({b,v,-1});
        }
    }
    events.push_back({stop,0,0});
    std::sort(events.begin(),events.end(),[](auto l,auto r){return l.x<r.x;});
    std::map<std::pair<int,int>,List> by_key;
    Q previous{0,1};
    for(std::size_t i=0;i<events.size();){
        const auto x=events[i].x;
        if(active.empty()) throw std::runtime_error("Dirichlet coverage failed");
        if(active.size()<=2){
            const int first=*active.begin(), last=active.size()==1?0:*active.rbegin();
            append(by_key[{first,last}],{previous,x});
        }
        while(i<events.size() && equal(events[i].x,x)){
            const auto e=events[i++];
            if(e.change==1) active.insert(e.v);
            if(e.change==-1) active.erase(e.v);
        }
        previous=x;
    }
    return by_key;
}

// The three types of cells that become free after deleting r and s.
List free_after_pair(std::map<std::pair<int,int>,List>& buckets,int r,int s){
    List raw=buckets[{r,0}];
    for (auto key : {std::make_pair(s,0),std::make_pair(r,s)}) {
        const auto& extra=buckets[key];
        raw.insert(raw.end(),extra.begin(),extra.end());
    }
    std::sort(raw.begin(),raw.end(),[](auto l,auto t){return l.a<t.a;});
    List result;
    for(auto cell:raw) append(result,cell);
    // Earliest failures tend to occur on the longest interval.
    std::sort(result.begin(),result.end(),[](auto l,auto t){return length(t)<length(l);});
    return result;
}

bool one_band_covers(Interval i,Int w,int n){
    // Need w*b-1 <= n*j <= w*a+1. The largest permissible j suffices.
    Int j=(w*i.a.p+i.a.d)/(Int(n)*i.a.d);
    return w*i.b.p <= (Int(n)*j+1)*i.b.d;
}

// Subtract the closed bands of w from an open interval. Only the open
// residual pieces matter: a strict witness always has a neighbourhood.
void subtract(Interval i,Int w,int n,List& remaining){
    Int j=(w*i.a.p)/(Int(n)*i.a.d);
    Q cursor=i.a;
    while(cursor<i.b){
        Q a{Int(n)*j-1,w}, b{Int(n)*j+1,w};
        if(!(a<i.b)) break;
        if(cursor<a) remaining.push_back({cursor,a});
        if(cursor<b) cursor=b;
        ++j;
    }
    if(cursor<i.b) remaining.push_back({cursor,i.b});
}

void run(int n){
    auto buckets=cells(n);
    int hits=0;
    Int first_count=0,second_count=0,largest_first=0,largest_second=0;
    for(int r=1;r<n;++r) for(int s=r+1;s<n;++s){
        const List free=free_after_pair(buckets,r,s);
        if(free.empty()) throw std::runtime_error("Empty deletion region");
        const Q longest=length(free.front());
        const Int first_cap=4*n*longest.d/((n-4)*longest.p);
        largest_first=std::max(largest_first,first_cap);
        for(Int u=n;u<=first_cap;++u){
            ++first_count;
            Int cap=std::numeric_limits<Int>::max();
            List residual;
            for(const auto& cell:free){
                const auto before=residual.size();
                subtract(cell,u,n,residual);
                for(std::size_t i=before;i<residual.size();++i){
                    const Q gap=length(residual[i]);
                    cap=std::min(cap,2*gap.d/gap.p);
                }
                if(cap<=u) break;
            }
            if(cap<=u) continue;
            if(residual.empty()) throw std::runtime_error("Single repair covered two holes");
            // Sort by lengths only before the first subtraction (in free_after_pair).
            // Post-subtraction lengths can have larger denominators.
            largest_second=std::max(largest_second,cap);
            for(Int w=u+1;w<=cap;++w){
                ++second_count;
                if(std::all_of(residual.begin(),residual.end(),[=](auto i){return one_band_covers(i,w,n);})){
                    std::cout<<"hit,"<<n<<','<<r<<','<<s<<','<<u<<','<<w<<'\n';
                    ++hits;
                }
            }
        }
    }
    std::cout<<"count,"<<n<<','<<first_count<<','<<second_count<<','<<largest_first
             <<','<<largest_second<<','<<hits<<'\n';
}
int main(int argc,char**argv){
    if(argc!=3){std::cerr<<"Usage: census FIRST_N LAST_N\n";return 2;}
    const int a=std::stoi(argv[1]),b=std::stoi(argv[2]);
    if(a<8 || b>300 || b<a){std::cerr<<"Expected 8 <= FIRST_N <= LAST_N <= 300\n";return 2;}
    for(int n=a;n<=b;++n) run(n);
}
