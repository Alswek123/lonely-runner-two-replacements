# Equality cases in the Lonely Runner problem

The latest proof draft classifies the equality cases for **six and seven distinct positive integer speeds of arbitrary size**. With six moving runners and a stationary reference, the unique primitive set with maximum loneliness $1/7$ is $\{1,2,3,4,5,6\}$. With seven moving runners there are exactly three at $1/8$. Every nonprimitive example is a common integer dilation. The individual examples were already known. The draft gives a computer-assisted argument for completeness and awaits external mathematical review. [Proof and reproduction instructions](equality/) · [PDF](equality/equality.pdf) · [中文说明](README.zh-CN.md)

This also closes the two remaining small sizes in the four-replacement study. For the complete four-replacement result, use the new corollary together with the earlier large-size argument. [Four-replacement proof and certificate](four-replacements/)

## Earlier work on two replacement speeds

[Read the proof draft (PDF)](runner.pdf) · [LaTeX source](runner.tex) · [中文说明](README.zh-CN.md)

Replace exactly two entries of the integer speed set $\{1,\ldots,n-1\}$ with two new integers $n\le u<w$. Call the resulting set *tight* if the greatest distance one can simultaneously put between every runner and the origin on the unit circle is exactly $1/n$.

The draft gives the following bounds for tight sets as $n$ grows.

$$u=O(n\log\log n),\qquad w=O(n\log n).$$

The first bound has the optimal order along an infinite sequence of tight examples. For the second bound, Section 5 uses several rational times associated with the same deleted speed to force divisibility. That draft left the optimal growth rate for $w$ open.

The question about the maximum replacement speed appears in Yuhan Zhang's Question 8.2. This note treats the case of exactly two replacements. The construction for the lower bound specializes the simultaneous accelerations of Goddyn and Wong. [Zhang, Section 8](https://arxiv.org/html/2608.13599v3#S8) · [Goddyn and Wong, *Tight Instances of the Lonely Runner*](https://zenodo.org/records/8275490)

## Review status

This research draft and its software were prepared using OpenAI Codex. The mathematical proof awaits independent review. Corrections, counterexamples and references to earlier equivalent results are welcome in [Issues](https://github.com/Alswek123/lonely-runner-two-replacements/issues).

## Reproduce the computation

The census covers every pair of deletions for $8\le n\le 300$. Each candidate's speed limit follows from a rational interval bound. Over 4,455,065 deletion pairs, three tight configurations survive, one each at $n=8,74,284$. All three examples were previously known. [Surviving speed sets](tight_sets.csv) · [Counts by size](census_by_n.csv)

Run this short check from the repository root with Python 3. It independently enumerates every size from 8 through 18 using exact fractions, checks each of the three survivors directly, and tests three parameter choices for the lower-bound construction.

```sh
python crosscheck.py
```

To compile the C++17 census and compare its output with both archived CSV tables, run

```sh
g++ -std=c++17 -O3 census.cpp -o census
python replay.py ./census 8 300
```

For a short replay containing a known example, replace `8 300` with `70 75`. On Windows, pass `.\census.exe` as the executable path.
