// Independent finite filter. Exact joint mass on a common interval in x.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
using Z=std::int64_t;
using Row=std::vector<int>;
int inv(int a,int m){
    int r=m,b=a,x=0,y=1;
    while(b){int q=r/b,t=r-q*b;r=b;b=t;t=x-q*y;x=y;y=t;}
    return (x%m+m)%m;
}
bool possible(int n,int r,int lo,int hi,int den,bool forbid){
    int s=n-r;
    if(2*r<n && 2*lo<den)return true;
    int stop=std::min(n-1,int(Z(s)*den/hi));
    Row fs,ps;
    for(int f=s;f<=stop;++f)if(std::gcd(f,r)==1){
        fs.push_back(f);ps.push_back(r-inv(f%r,r));
        if(fs.size()==64)break;
    }
    if(fs.empty())return true;
    Z cap=8*Z(n)*r*den/((n-8)*(hi-lo));
    Z maxj=(Z(r)*den+cap*hi)/(Z(n)*den);
    if(cap>3000000 || maxj>100)return true;
    Row contacts(r);
    for(int f:fs)for(int j=1;j<=maxj;++j)++contacts[Z(j)*f%r];
    Z period=Z(n)*r*den,radius=Z(r)*den;
    for(int rem=0;rem<r;++rem){
        if(rem && 4*contacts[rem]<int(fs.size()))continue;
        Z v=n+(rem-(n%r)+r)%r;
        for(;v<=cap;v+=r){
            if(!rem && forbid && 2*(v/r)*s<=n)continue;
            Z width=v*(hi-lo),total=0;
            for(int p:ps){
                Z start=Z(n)*den*(v*p%r)+v*lo, end=start+width;
                for(Z j=(start+radius)/period;j*period-radius<end;++j){
                    Z l=std::max(start,j*period-radius),h=std::min(end,j*period+radius);
                    if(h>l)total+=h-l;
                }
            }
            if(4*total>=Z(ps.size())*width)return true;
        }
    }
    return false;
}
int main(int argc,char**argv){
    if(argc!=2)return 2;
    std::ifstream in(argv[1]);if(!in)return 2;
    int n,a,b;
    while(in>>n>>a>>b){
        if(n<9||n>4000)return 2;
        Row good(a),bad;
        for(int& r:good)if(!(in>>r))return 2;
        for(int i=0;i<b;++i){
            int r;if(!(in>>r))return 2;
            if(n<75||r==1){bad.push_back(r);continue;}
            int s=n-r;
            if(!possible(n,r,1,2,3,false)||!possible(n,r,1,3,4,false)||
               !possible(n,r,1,2,2,false)||!possible(n,r,3,4,6,false)||
               !possible(n,r,1,2,4,false))continue;
            int fmax=n/4;
            if(2*r>=n && s<=fmax && !possible(n,r,s,2*s,2*fmax,true)){
                good.push_back(r);continue;
            }
            bad.push_back(r);
        }
        if(good.size()+bad.size()<4 || bad.empty())continue;
        std::cout<<n<<' '<<good.size()<<' '<<bad.size();
        for(int r:good)std::cout<<' '<<r;
        for(int r:bad)std::cout<<' '<<r;
        std::cout<<'\n';
    }
}
