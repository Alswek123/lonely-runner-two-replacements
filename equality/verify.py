"""Replay the prime certificate and the exact continuous-time equality checks."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from math import gcd, isqrt, prod
from pathlib import Path
import subprocess

MODELS = {
    7: [(1, 2, 3, 4, 5, 6)],
    8: [(1, 2, 3, 4, 5, 6, 7), (1, 2, 3, 4, 5, 7, 12), (1, 4, 5, 6, 7, 11, 13)],
}
PRIMES = {
    7: (373, 379, 383, 389, 397, 401, 409, 419, 433, 439, 443, 449, 463, 467, 479, 491),
    8: (367, 379, 383, 389, 397, 419, 431, 463, 487, 499, 503, 509, 521, 541, 547,
        557, 563, 569, 587, 593, 599, 601, 607),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fold(x, p):
    r = x % p
    return min(r, p-r)


def exact_lr(model):
    # Each change from a positive to a negative slope of the lower envelope
    # occurs when (v+w)t is integral, including v=w at a triangle peak.
    times = {Fraction(j, v+w) for v, w in combinations_with_replacement(model, 2)
             for j in range(v+w)}
    return max(min(Fraction(fold(t.numerator*v, t.denominator), t.denominator)
                   for v in model) for t in times)


def reference_rows(n, p):
    return {tuple(sorted(fold(v*pow(b, -1, p), p) for v in model))
            for model in MODELS[n] for b in model}


def parse_rows(text):
    rows = [tuple(map(int, line.split())) for line in text.splitlines()]
    require(len(rows) == len(set(rows)), 'duplicate row')
    return set(rows)


def test(exe, n, p):
    raw = subprocess.check_output([exe, str(n), str(p)], text=True)
    actual = parse_rows(raw)
    expected = reference_rows(n, p)
    require(actual == expected,
            f'n={n}, p={p}, missing={expected-actual}, extra={actual-expected}')
    for row in actual:
        require(all(min(fold(t*v,p) for v in row)*n <= p for t in range(1,p)),
                f'non-cover {row} at {p}')
    return n, p, len(actual)


def controls(exe):
    for n in (7, 8):
        p = 23
        all_times = set(range(1,p))
        classes = list(range(1,p//2+1))
        bands = {v: {t for t in all_times if n*fold(t*v,p) <= p} for v in classes}
        minimal = set()
        for count in range(n-1):
            for extra in combinations(classes[1:], count):
                row = frozenset((1, *extra))
                if any(previous <= row for previous in minimal):
                    continue
                if set.union(*(bands[v] for v in row)) == all_times:
                    minimal.add(row)
        reported = {frozenset(row) for row in parse_rows(subprocess.check_output(
            [exe, str(n), str(p)], text=True))}
        require(minimal <= reported and all(any(m <= r for m in minimal) for r in reported),
                f'brute-force control at n={n}')
    n, p = 7, 367
    rows = parse_rows(subprocess.check_output([exe, str(n), str(p)], text=True))
    extras = rows - reference_rows(n, p)
    require(bool(extras) and all(all(min(fold(t*v,p) for v in row)*n <= p
                                     for t in range(1,p)) for row in extras),
            'unlisted cover expected for the negative control')
    print('PASS: direct subset checks at p=23 and unexpected covers at n=7,p=367.',
          flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path, help='compiled covers or crosscheck')
    parser.add_argument('--jobs', type=int, default=1)
    args = parser.parse_args()
    exe = args.executable.resolve(strict=True)
    require(args.jobs >= 1, 'jobs must be positive')
    controls(exe)
    for n, primes in PRIMES.items():
        k = n*(n-1)//2
        size_bound = k**(n-2)
        require(len(set(primes)) == len(primes), 'repeated prime')
        for p in primes:
            require(p >= 2 and all(p%d for d in range(2,isqrt(p)+1)), f'nonprime {p}')
            require(all(k*max(model) < p for model in MODELS[n]), 'lifting threshold')
        require(prod(primes) > size_bound**(n-1), 'insufficient prime product')
        print(f'n={n}: W <= {size_bound}, product = {prod(primes)} > {size_bound**(n-1)}',
              flush=True)
        for model in MODELS[n]:
            value = exact_lr(model)
            require(gcd(*model) == 1 and value == Fraction(1,n), f'model {model}: {value}')
            print(f'  {model}: LR = {value}')
        four = []
        for model in MODELS[n]:
            for d in range(1,n):  # At least one retained speed lies in [n-1].
                dilated = tuple(d*x for x in model)
                if sum(v >= n for v in dilated) == 4:
                    four.append((d, model))
        expected = [(3, MODELS[7][0])] if n == 7 else [(2, m) for m in MODELS[8][:2]]
        require(four == expected, 'four-replacement corollary')
    cases = [(n, p) for n, primes in PRIMES.items() for p in primes]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(test, exe, n, p) for n, p in cases]
        for future in futures:
            n, p, count = future.result()
            print(f'{n} {p} {count}', flush=True)
    print(f'PASS: {len(cases)} prime certificates, exact template values and four-replacement cases.')


if __name__ == '__main__':
    main()
