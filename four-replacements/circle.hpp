// Exact scaled-time cell decomposition shared by the final four-speed enumeration.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using Z=std::int64_t;
using Wide=__int128;
struct Q { Z p,d; }; // scaled time x=n*t
struct Cell { Q l,r; };
struct Length { Wide p,d; };
bool less(Q x,Q y){return Wide(x.p)*y.d<Wide(y.p)*x.d;}
bool eq(Q x,Q y){return !less(x,y)&&!less(y,x);}
Length length(Cell c){return {Wide(c.r.p)*c.l.d-Wide(c.l.p)*c.r.d,Wide(c.l.d)*c.r.d};}
bool longer(Length x,Length y){return x.p*y.d>y.p*x.d;}
using List=std::vector<Cell>;
using Key=std::array<int,4>;
using Buckets=std::map<Key,List>;
struct Event{Q x;int v;bool add;};
void append(List& v,Cell c){
    if(!less(c.l,c.r))return;
    if(!v.empty()&&eq(v.back().r,c.l))v.back().r=c.r;
    else v.push_back(c);
}
Buckets buckets(int n){
    std::vector<Event> events;
    std::set<int> active;
    Q stop{n,2};
    for(int v=1;v<n;++v){
        active.insert(v);
        for(int j=0;j<=v/2+1;++j){
            Q l{Z(n)*j-1,v},r{Z(n)*j+1,v};
            if(l.p>0&&less(l,stop))events.push_back({l,v,true});
            if(less(r,stop))events.push_back({r,v,false});
        }
    }
    events.push_back({stop,0,false});
    std::sort(events.begin(),events.end(),[](Event a,Event b){return less(a.x,b.x);});
    Q prev{0,1};
    Buckets by;
    for(std::size_t j=0;j<events.size();){
        Q x=events[j].x;
        if(active.empty())throw std::runtime_error("baseline gap");
        if(active.size()<=4){
            Key key{}; std::copy(active.begin(),active.end(),key.begin());
            append(by[key],{prev,x});
        }
        while(j<events.size()&&eq(events[j].x,x)){
            auto ev=events[j++];
            if(ev.add)active.insert(ev.v);else active.erase(ev.v);
        }
        prev=x;
    }
    return by;
}
List uncovered(const Buckets& by, Key removed){
    List raw,merged;
    for(int mask=1;mask<16;++mask){
        Key key{};int used=0;
        for(int j=0;j<4;++j)if(mask&(1<<j))key[used++]=removed[j];
        auto found=by.find(key);
        if(found!=by.end())raw.insert(raw.end(),found->second.begin(),found->second.end());
    }
    std::sort(raw.begin(),raw.end(),[](Cell a,Cell b){return less(a.l,b.l);});
    for(auto c:raw)append(merged,c);
    return merged;
}
bool covers(Cell c,Z v,int n){
    Z j=(v*c.l.p+c.l.d)/(n*c.l.d);
    return Wide(v)*c.r.p<=Wide(n*j+1)*c.r.d;
}
void subtract(Cell c,Z v,int n,List& out){
    Z j=(v*c.l.p)/(n*c.l.d);
    Q pos=c.l;
    while(less(pos,c.r)){
        Q l{n*j-1,v},r{n*j+1,v};
        if(!less(l,c.r))break;
        if(less(pos,l))out.push_back({pos,l});
        if(less(pos,r))pos=r;
        ++j;
    }
    if(less(pos,c.r))out.push_back({pos,c.r});
}
