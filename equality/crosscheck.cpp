// A second enumerator. It stores both incidences of the time/speed relation,
// uses value-owned candidate sets, and solves the last two slots by incidence.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <vector>

// All certificate primes are below 640. Bits are indexed from 1.
struct Bits {
    std::array<std::uint64_t, 5> word{};
    void put(int i) { word[i/64] |= std::uint64_t(1) << (i%64); }
    void erase(int i) { word[i/64] &= ~(std::uint64_t(1) << (i%64)); }
    int size() const {
        int total = 0;
        for (auto w : word) total += __builtin_popcountll(w);
        return total;
    }
    template<class F> void each(F f) const {
        for (int i = 0; i < 5; ++i) {
            auto w = word[i];
            while (w) {
                f(64*i + __builtin_ctzll(w));
                w &= w-1;
            }
        }
    }
    Bits minus(const Bits& b) const {
        Bits result = *this;
        for (int i=0; i<5; ++i) result.word[i] &= ~b.word[i];
        return result;
    }
    Bits intersect(const Bits& b) const {
        Bits result = *this;
        for (int i=0; i<5; ++i) result.word[i] &= b.word[i];
        return result;
    }
};

struct Enumerator {
    int n, p, m;
    std::vector<Bits> times, speeds;
    std::set<std::vector<int>> output;
    std::vector<int> prefix{1};

    Enumerator(int runners,int prime) : n(runners),p(prime),m(p/2),
        times(m+1),speeds(m+1) {
        // Generate bad times by inverting v on the residues [-floor(p/n),floor(p/n)].
        for (int v=1; v<=m; ++v) {
            int inv=1;
            while (v*inv%p != 1) ++inv;
            for (int b=1; b<=p/n; ++b) {
                int t=b*inv%p;
                t=std::min(t,p-t);
                times[v].put(t);
                speeds[t].put(v);
            }
        }
    }

    void record() {
        auto row=prefix;
        std::sort(row.begin(),row.end());
        output.insert(row);
    }

    void search(Bits gaps, Bits candidates, int left) {
        int missing=gaps.size();
        if (!missing) { record(); return; }
        if (left==2) {
            // At least one of the pair covers the selected gap.
            int pivot=0;
            gaps.each([&](int t){ pivot=t; });
            candidates.intersect(speeds[pivot]).each([&](int a){
                Bits residual=gaps.minus(times[a]);
                prefix.push_back(a);
                if (!residual.size()) record();
                else {
                    Bits possible=candidates;
                    residual.each([&](int t){ possible=possible.intersect(speeds[t]); });
                    possible.each([&](int b){ prefix.push_back(b); record(); prefix.pop_back(); });
                }
                prefix.pop_back();
            });
            return;
        }
        // Rank by the number of gaps hit. If the sum of the next 'left'
        // weights is too small, no suffix of this ordering contains a cover.
        std::vector<std::pair<int,int>> order;
        candidates.each([&](int v){
            int hits=times[v].intersect(gaps).size();
            if (hits) order.emplace_back(hits,v);
        });
        std::sort(order.rbegin(),order.rend());
        std::vector<int> sums(order.size()+1);
        for (std::size_t j=0;j<order.size();++j) sums[j+1]=sums[j]+order[j].first;
        for (std::size_t j=0;j<order.size();++j) {
            auto end=std::min(order.size(),j+left);
            if (sums[end]-sums[j]<missing) break;
            int v=order[j].second;
            candidates.erase(v);
            prefix.push_back(v);
            search(gaps.minus(times[v]),candidates,left-1);
            prefix.pop_back();
        }
    }
    void run() {
        Bits full;
        for (int i=1;i<=m;++i) full.put(i);
        Bits candidates=full;
        candidates.erase(1);
        search(full.minus(times[1]),candidates,n-2);
        for (const auto& row:output) {
            for (int value:row) std::cout << value << ' ';
            std::cout << '\n';
        }
    }
};

int main(int argc,char** argv) {
    if (argc!=3) return 2;
    int n=std::stoi(argv[1]),p=std::stoi(argv[2]);
    if ((n!=7 && n!=8) || p<3 || p>=640 || !(p%2)) return 2;
    for (int d=3;d*d<=p;d+=2) if (p%d==0) return 2;
    Enumerator(n,p).run();
    return std::cout ? 0 : 1;
}
