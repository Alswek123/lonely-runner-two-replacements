"""Independent exact checks with Python's Fraction. No third-party dependencies."""
from fractions import Fraction as F
from itertools import combinations
from functools import lru_cache
from math import floor, gcd, isqrt, prod


@lru_cache(maxsize=1024)
def good(n, v):
    """Good time intervals in ordinary t-coordinates on the half circle."""
    result = []
    for k in range(v // 2 + 1):
        a = F(n * k + 1, n * v)
        b = min(F(1, 2), F(n * (k + 1) - 1, n * v))
        if a < b:
            result.append((a, b))
    return result


def overlap(left, right):
    output = []
    i = j = 0
    while i < len(left) and j < len(right):
        a, b = left[i]
        c, d = right[j]
        if max(a, c) < min(b, d):
            output.append((max(a, c), min(b, d)))
        if b <= d:
            i += 1
        if d <= b:
            j += 1
    return output


def eligible(n, speeds):
    windows = [(F(0), F(1, 2))]
    for v in speeds:
        windows = overlap(windows, good(n, v))
    return windows


def exact_census(n):
    hits = []
    base = set(range(1, n))
    for r, s in combinations(sorted(base), 2):
        windows = eligible(n, sorted(base - {r, s}))
        length = max(b - a for a, b in windows)
        cap = floor(4 / ((n - 4) * length))
        for u in range(n, cap + 1):
            remaining = overlap(windows, good(n, u))
            if not remaining:
                raise ValueError("A single inserted speed covered two deleted speeds")
            longest = max(b - a for a, b in remaining)
            for w in range(u + 1, floor(2 / (n * longest)) + 1):
                if not overlap(remaining, good(n, w)):
                    hits.append((n, r, s, u, w))
    return hits


def check_gw_example(n, r, s, u, w):
    # Test the original inequalities independently, including equality at t=1/n.
    speeds = sorted((set(range(1, n)) - {r, s}) | {u, w})
    assert not eligible(n, speeds)
    assert min(min(v % n, n - v % n) for v in speeds) == 1
    if n != 8:
        for v, new in ((r, u), (s, w)):
            m, rem = divmod(new, v)
            assert rem == 0 and m >= 2
            assert all(gcd(v, b) > 1 for b in range(n - v, m * (n - v)))


def check_lower_construction(m):
    def primorial(limit):
        return prod(p for p in range(2, limit + 1)
                    if all(p % d for d in range(2, isqrt(p) + 1)))

    s = primorial(2 * m) + 2
    q = primorial(m * s)
    n = q + s
    assert n <= m * q < m * (n - 2)
    for d in (2, s):
        assert all(gcd(n - d, b) > 1 for b in range(d, m * d))
    return n, m * q, m * (n - 2)


if __name__ == "__main__":
    expected = [(8, 2, 3, 11, 13), (74, 70, 72, 140, 144),
                (284, 280, 282, 560, 564)]
    hits = [row for n in range(8, 19) for row in exact_census(n)]
    assert hits == expected[:1], hits
    for row in expected:
        check_gw_example(*row)
    examples = [check_lower_construction(m) for m in (2, 3, 4)]
    assert examples[0] == (30038, 60060, 60072)
    print("Independent rational census for n=8..18 and all three tight witnesses passed.")
    print("Lower-construction divisibility checks passed for m=2, 3, 4.")
