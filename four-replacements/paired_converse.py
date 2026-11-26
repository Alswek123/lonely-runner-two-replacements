"""Construct the blocking integers in the converse proof. Standard-library Python."""
from functools import lru_cache
from math import gcd, isqrt, prod
import csv

PRIME_LADDER = (5,17,53,157,467,1399,4177,12527,37579,112691,338033,1014089)


@lru_cache(None)
def primes_of(n):
    rest, ans = n, []
    p = 2
    while p*p <= rest:
        if rest%p == 0:
            ans.append(p)
            while rest%p == 0:
                rest //= p
        p += 1 if p==2 else 2
    if rest>1:
        ans.append(rest)
    return tuple(ans)


def chain_data(e):
    p = next(z for z in range(2,e+2) if primes_of(z)==(z,) and e%z)
    chain = [1]
    for z in range(2,e+2):
        if gcd(z,e)!=1:
            continue
        if z>2*chain[-1]:
            chain = []
        chain.append(z)
        if p*chain[0] <= 2*z:
            return p,chain
    return None


def check_finite_lemmas():
    # Each chain repeats after multiplication by p, with all ratios at most two.
    cases = 0
    for e in range(4,5018,2):
        if e==6:
            continue
        p,chain = chain_data(e)
        assert primes_of(p)==(p,) and gcd(p,e)==1
        assert all(v<=e+1 and gcd(v,e)==1 for v in chain)
        assert all(y<=2*x for x,y in zip(chain,chain[1:]))
        assert p*chain[0]<=2*chain[-1]
        cases += 1
    assert 2*13*(6*2**5+1) == 5018
    assert 2*17*(7*2**6+1) < prod((2,3,5,7,11,13))
    end = 1
    for q in PRIME_LADDER:
        assert all(q%d for d in range(2,isqrt(q)+1))
        t=1
        while 2*t<=q:
            t*=3
        assert 3*t<2*q<4*t and q//3+1<=end+1<=q
        end=q
    assert end>=396739
    return cases


def third_block(s,a):
    # The archived rows have a inside the printed certificate.
    q = next(q for q in PRIME_LADDER if a<=q<3*a)
    t = 1
    while 2*t<=q:
        t*=3
    for base in (t,q):
        y=base
        while y<s:
            y*=3
        if y<2*s:
            assert y%3 or y%9==0 or 7*y<9*s
            return y
    raise AssertionError((s,a))


def check_archive(file):
    counts = {3:0, 5:0}
    with open(file,newline='') as stream:
        for s,f,m,a,h,old_witness in (map(int,row) for row in csv.reader(stream)):
            d=s-a
            g=gcd(d,f)
            e=d//g
            assert e%2==0 and 2<=e<f//g<=h and g>=h*a
            assert all(g%p for p in range(2,h+1) if primes_of(p)==(p,))
            assert a*h<=f and h*d%f==0
            if h==3:
                y=third_block(s,a)
                counts[3]+=1
            else:
                assert h>=5
                y=next(y for y in range(s,2*s) if gcd(y,e)==1 and max(primes_of(y))<=h
                       and (y%3!=0 or y%9==0))
                counts[5]+=1
            assert s<=y<2*s and y>h and gcd(y,d)==1
    assert sum(counts.values())==1337
    return counts


def gcd_tests(n,shape):
    return all(gcd(f,n-s)>1 for s,m in shape for f in range(s,m*s))


def certificate(n,shape):
    """Return None for a GW pattern, else exact t=numerator/denominator with all distances >1/n.
    Assumes distinct valid replacements, and 2*max(m*s)<=n."""
    shape=tuple(sorted(shape))
    if not shape or len({s for s,m in shape})!=len(shape):
        raise ValueError('Supply distinct deleted offsets.')
    if any(not (1<=s<n and m>=2 and m*(n-s)>=n) for s,m in shape):
        raise ValueError('Invalid replacement.')
    speeds=[m*(n-s) for s,m in shape]
    if len(set(speeds))!=len(speeds) or 2*max(m*s for s,m in shape)>n:
        raise ValueError('This test requires distinct speeds and 2*max(m*s)<=n.')
    for i,(s,m) in enumerate(shape):
        r=n-s
        units=[f for f in range(s,m*s) if gcd(f,r)==1]
        if not units:
            continue
        for f in units:
            if any(h*a<=f and h*(s-a)%f==0 for a,h in shape[:i]):
                continue
            p=(-pow(f,-1,r))%r
            ell=(f*p+1)//r
            w=max(speeds)
            num=2*w*(n*ell-1)-1
            den=2*n*f*w
            div=gcd(num,den)
            return num//div,den//div
        raise AssertionError('Every exceptional unit has an arithmetic source.')
    return None


def check_time(n,shape,t):
    a,b=t
    deleted={n-s for s,h in shape}
    speeds=set(range(1,n)) - deleted | {h*(n-s) for s,h in shape}
    return len(speeds)==n-1 and all(n*min(v*a%b,(-v*a)%b)>b for v in speeds)
