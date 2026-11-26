// List every inclusion-minimal cover reachable with at most n-1 classes,
// or a shorter cover if one exists. Output every full cover when no shorter
// cover exists. Class 1 is fixed by the multiplicative action modulo p.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using Mask = std::vector<std::uint64_t>;

class Search {
    int n, p, m, words;
    std::vector<Mask> band;
    std::vector<bool> allowed;
    std::vector<int> chosen{1};

    void visit(const Mask& gaps, int count, int slots) {
        if (count == 0) {
            auto row = chosen;
            std::sort(row.begin(), row.end());
            for (int v : row) std::cout << v << ' ';
            std::cout << '\n';
            return;
        }
        if (!slots) return;
        Mask next(words);
        std::vector<int> tried;
        for (int v = 2; v <= m; ++v) {
            if (!allowed[v]) continue;
            int hits = 0;
            for (int j = 0; j < words; ++j)
                hits += __builtin_popcountll(gaps[j] & band[v][j]);
            if (hits * slots < count) continue;
            for (int j = 0; j < words; ++j) next[j] = gaps[j] & ~band[v][j];
            allowed[v] = false;
            chosen.push_back(v);
            visit(next, count - hits, slots - 1);
            chosen.pop_back();
            tried.push_back(v);
            // Leave v excluded for the remaining siblings.
        }
        for (int v : tried) allowed[v] = true;
    }

public:
    Search(int runners, int prime)
        : n(runners), p(prime), m(p/2), words((m+63)/64),
          band(m+1, Mask(words)), allowed(m+1, true) {}

    void run() {
        Mask gaps(words);
        int count = 0;
        for (int t = 1; t <= m; ++t) {
            int pos = (t-1)/64;
            auto bit = std::uint64_t(1) << ((t-1)%64);
            for (int v = 1; v <= m; ++v) {
                int r = t*v % p;
                if (n * std::min(r, p-r) <= p) band[v][pos] |= bit;
            }
            if (!(band[1][pos] & bit)) {
                gaps[pos] |= bit;
                ++count;
            }
        }
        visit(gaps, count, n-2);
    }
};

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    int n = std::stoi(argv[1]), p = std::stoi(argv[2]);
    if ((n != 7 && n != 8) || p < 3 || p > 4093 || p % 2 == 0) return 2;
    for (int d = 3; d*d <= p; d += 2) if (p % d == 0) return 2;
    Search(n,p).run();
    return std::cout ? 0 : 1;
}
