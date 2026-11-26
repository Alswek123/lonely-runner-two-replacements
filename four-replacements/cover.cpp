#include "circle.hpp"
#include <numeric>
using Hit=std::array<Z,4>;
using Hits4=std::set<Hit>;
struct Stats{Z nodes=0;std::array<Z,4> maxima{};};

bool enough(Cell c,int n,Z v,int left){
    const Z den=c.l.d*c.r.d,lo=v*c.l.p*c.r.d,hi=v*c.r.p*c.l.d;
    Z mass=0;
    for(Z j=(lo+den)/(n*den);(n*j-1)*den<hi;++j){
        Z a=std::max(lo,(n*j-1)*den),b=std::min(hi,(n*j+1)*den);
        if(b>a)mass+=b-a;
    }
    return left*mass>=hi-lo;
}
void finish(int n,const List& holes,const std::vector<int>& centres,std::vector<Z>& vs,
            Hits4& hits,Stats& stats,int initial_bound=0){
    if(holes.empty())throw std::runtime_error("earlier insertions already cover everything");
    int left=4-int(vs.size());
    Cell probe=holes[0];Length len=length(probe);
    for(Cell c:holes)if(longer(length(c),len)){probe=c;len=length(c);}
    Wide upper;
    if(left==1)upper=2*len.d/len.p;
    else if(n>2*left)upper=2*left*n*len.d/((n-2*left)*len.p);
    else if(initial_bound>0 && vs.empty())upper=initial_bound;
    else throw std::runtime_error("first-speed bound needed");
    if(upper>20000)throw std::runtime_error("speed guard exceeded");
    Z cap=Z(upper),step=1;
    ++stats.nodes;
    stats.maxima[left-1]=std::max(stats.maxima[left-1],cap);
    if(left==1){
        for(int q:centres){
            bool already=false;for(Z v:vs)if(v%q==0)already=true;
            if(!already)step=std::lcm(step,Z(q));
            if(step>cap)return;
        }
    }
    for(Z v=((n+step-1)/step)*step;v<=cap;v+=step){
        if(n<=8 && !vs.empty() && v<vs.front())continue;
        if(std::find(vs.begin(),vs.end(),v)!=vs.end())continue;
        if(left==1){
            if(covers(probe,v,n) && std::all_of(holes.begin(),holes.end(),[&](Cell c){return covers(c,v,n);})){
                Hit hit{vs[0],vs[1],vs[2],v};
                std::sort(hit.begin(),hit.end());hits.insert(hit);
            }
        }else{
            if((n>2*left || !vs.empty()) && !enough(probe,n,v,left))continue;
            List more;
            for(Cell c:holes)subtract(c,v,n,more);
            vs.push_back(v);finish(n,more,centres,vs,hits,stats);vs.pop_back();
        }
    }
}
int main(int argc,char**argv){
    if(argc!=4)return 2;
    std::ifstream in(argv[1]);if(!in)return 2;
    int lo=std::stoi(argv[2]),hi=std::stoi(argv[3]),n,initial;
    std::map<int,std::map<Key,int>> cases;
    Key key;
    while(in>>n>>initial>>key[0]>>key[1]>>key[2]>>key[3]){
        if(n>=lo && n<=hi)cases[n][key]=initial;
    }
    if(!in.eof())return 2;
    for(const auto& [n,group]:cases){
        if(n>100 || n<7)throw std::logic_error("outside checked integer range");
        auto table=buckets(n);
        Stats stats;Z count=0;
        for(const auto& [removed,bound]:group){
            auto holes=uncovered(table,removed);
            std::vector<int> centres;
            for(int q=2;q<n;++q){
                bool all=true;
                for(int r=q;r<n;r+=q)if(std::find(removed.begin(),removed.end(),r)==removed.end())all=false;
                if(all)centres.push_back(q);
            }
            std::vector<Z> vs;Hits4 hits;
            finish(n,holes,centres,vs,hits,stats,bound);
            for(auto hit:hits){
                std::cout<<"H "<<n;for(int r:removed)std::cout<<' '<<r;
                for(Z v:hit)std::cout<<' '<<v;
                std::cout<<'\n';++count;
            }
        }
        std::cout<<"C "<<n<<' '<<group.size()<<' '<<stats.nodes;
        for(Z c:stats.maxima)std::cout<<' '<<c;
        std::cout<<' '<<count<<std::endl;
    }
}
