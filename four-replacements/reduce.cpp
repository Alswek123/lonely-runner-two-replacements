// Finite deletion filter using only the two averaged inequalities in the paper.
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>
using Row=std::vector<int>;
int main(int argc,char**argv){
    if(argc!=2)return 2;
    const int limit=std::atoi(argv[1]);
    if(limit<9 || limit>11999)return 2;
    std::vector<Row> yes(limit+1),other(limit+1);
    for(int n=9;n<=std::min(24,limit);++n)other[n].push_back(1);
    for(int r=2;r<limit;++r){
        Row divisors{1};
        int rem=r;
        auto add_prime=[&](int p){
            int size=divisors.size();
            for(int j=0;j<size;++j)divisors.push_back(-p*divisors[j]);
        };
        for(int p=2;p<=rem/p;++p)if(rem%p==0){
            add_prime(p);do{rem/=p;}while(rem%p==0);
        }
        if(rem>1)add_prime(rem);
        auto units=[&](int x){int count=0;for(int d:divisors)count+=x/d;return count;};
        int phi=units(r),width=(60+int(divisors.size()))*r/phi+1;
        int last=limit-r;
        if(r>=width)last=std::min(last,std::max(999-r,3*width-1));
        for(int s=std::max(1,9-r);s<=last;++s){
            int n=r+s;
            int before=units(s-1);
            if(n>36){
                int h=units(std::min(n-1,4*s/3))-before;
                if(h*(n-36)>56*r)continue;
            }
            bool forced=false;
            if(n>28 && 2*r>=n && s<=n/4){
                int h=units(n/4)-before;
                forced=h*(n-28)>40*n;
            }
            (forced?yes[n]:other[n]).push_back(r);
        }
    }
    for(int n=9;n<=limit;++n){
        if(other[n].empty()||other[n].size()+yes[n].size()<4)continue;
        std::cout<<n<<' '<<yes[n].size()<<' '<<other[n].size();
        for(int r:yes[n])std::cout<<' '<<r;
        for(int r:other[n])std::cout<<' '<<r;
        std::cout<<'\n';
    }
}
