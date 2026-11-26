# Equality for six and seven integer speeds

Divide the distinct positive speeds by their common gcd. With one stationary reference runner, the proof draft gives the following complete list at equality. Every unlisted primitive set has a strictly better time.

| Moving speeds | Maximum loneliness | All primitive speed sets attaining the value |
|---|---|---|
| 6 | $1/7$ | $\{1,2,3,4,5,6\}$ |
| 7 | $1/8$ | $\{1,2,3,4,5,6,7\}$, $\{1,2,3,4,5,7,12\}$, $\{1,4,5,6,7,11,13\}$ |

The individual examples appear in earlier work. Goddyn and Wong specifically suggested uniqueness for six speeds in 2006. The contribution claimed here is completeness at arbitrary speed, subject to external review. The paper credits earlier dimensional reduction and prime-divisibility methods. [Manuscript](equality.pdf) · [LaTeX](equality.tex) · [Goddyn and Wong](https://zenodo.org/records/8275490) · [Rosenfeld](https://arxiv.org/abs/2509.14111) · [中文说明](README.zh-CN.md)

## Reproduce both enumerations

Use Python 3.9+ and a C++17 compiler supporting GCC/Clang bit-count builtins. From this directory run

```sh
c++ -O3 -std=c++17 covers.cpp -o covers
c++ -O3 -std=c++17 crosscheck.cpp -o crosscheck
python verify.py ./covers --jobs 4
python verify.py ./crosscheck --jobs 4
```

On Windows, compile with Clang or Zig C++, add the executable suffix and pass the corresponding path to the driver. The Python driver uses only the standard library.

Each executable takes the total number of runners and one prime as arguments. It emits the actual normalized modular covers. Neither enumerator contains the templates. The driver checks every row against the independently generated template orbits and evaluates every survivor on the full time grid. It also compares against brute-force subsets at prime 23, checks a prime with extra covers, verifies the size and prime-product inequalities, and recomputes the continuous-time maxima using exact fractions. Successful replays produce the same transcript. [Replay output](replay.txt)

## Four replacements

With six moving speeds the only four-replacement example is $3\{1,2,3,4,5,6\}$. With seven moving speeds there are two, obtained by doubling $\{1,2,3,4,5,6,7\}$ or $\{1,2,3,4,5,7,12\}$. Here a replacement removes an entry of $[n-1]$ and inserts a distinct integer at least $n$. The companion study treats the other sizes. [Four-replacement study](../four-replacements/)
