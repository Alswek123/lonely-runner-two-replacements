"""Run with Python 3.9+ and assertions enabled. For full regeneration add --full."""
import argparse
from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd, isqrt, prod
import os
from pathlib import Path
import random
import subprocess
import tempfile

from paired_converse import check_finite_lemmas

ROOT=Path(__file__).resolve().parent
FIRST=148431434


def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))


def first_four():
    start=64
    while start<=FIRST//2:
        last=min(5*start//4,FIRST//2)
        product=1
        for p in range(last,2*start):
            if prime(p):
                product*=p
                if product>FIRST:
                    break
        assert product>FIRST
        start=last+1
    candidates=[]
    for s in range(2,64):
        bases=[]
        for p in range(2,2*s):
            if not prime(p):
                continue
            value=p
            while value<s:
                value*=p
            if value<2*s:
                bases.append(p)
        if s+prod(bases)<=FIRST:
            candidates.append((s,set(bases)))
    # Assemble residues for each prime only after the four offsets are fixed.
    hits=[]
    for group in combinations(candidates,4):
        mods={}
        if any(len({s%p for s,ps in group if p in ps})>1
               for p in set.union(*(ps for s,ps in group))):
            continue
        for s,ps in group:
            mods.update({p:s%p for p in ps})
        modulus=prod(mods)
        x=sum(b*(modulus//p)*pow(modulus//p,-1,p) for p,b in mods.items())%modulus
        if x<=group[-1][0]:
            x+=((group[-1][0]-x)//modulus+1)*modulus
        if x<=FIRST:
            hits.append((x,tuple(s for s,ps in group)))
    assert hits==[(FIRST,(2,4,8,14))]
    first_units=[next(f for f in range(s,3*s) if gcd(f,FIRST-s)==1) for s in hits[0][1]]
    assert first_units==[5,9,17,29]
    return hits[0]


def small_intervals():
    for n,expected,ratio_cap in [(7,F(2,21),123),(8,F(1,14),23)]:
        lengths=[]
        for retained in combinations(range(1,n),n-5):
            walls={F(0),F(1)}
            for v in retained:
                for j in range(v+1):
                    walls.update(x for d in (-1,1) if 0<(x:=F(n*j+d,n*v))<1)
            ordered=sorted(walls)
            usable=[]
            for l,r in zip(ordered,ordered[1:]):
                mid=(l+r)/2
                if all(min(v*mid%1,-v*mid%1)>F(1,n) for v in retained):
                    usable.append(r-l)
            lengths.append(max(usable))
        assert min(lengths)==expected
        a=[F(1)]
        for j in (1,2,3):
            u=4-j
            c=F(2) if u==1 else F(2*u*n,n-2*u)
            den=n-2*j-2*j/(257*expected)
            assert den>0
            a.append(c*(sum(a)+(2*j+1)/(257*expected))/den)
        assert a[-1]<ratio_cap


def band_mass(n,v,left,right):
    """Direct integration in ordinary time, independent of the scaled interval code."""
    lo=int(left*v)-1
    hi=int(right*v)+2
    mass=F(0)
    for j in range(lo,hi+1):
        a=max(left,F(n*j-1,n*v))
        b=min(right,F(n*j+1,n*v))
        mass+=max(F(0),b-a)
    return mass/(right-left)


def average_checks():
    rng=random.Random(148431434)
    total=0
    for _ in range(900):
        n=rng.randrange(9,97)
        r=rng.randrange(2,n)
        s=n-r
        v=rng.choice([rng.randrange(n,6*n*r),r*rng.randrange(1+n//r,8*n)])
        if rng.randrange(2) and s<=n//4:
            stop=n//4
            b=F(s,stop)
            a=b/2
        else:
            stop=min(n-1,4*s//3)
            a,b=F(1,2),F(3,4)
        fs=[f for f in range(s,stop+1) if gcd(f,r)==1]
        if not fs:
            continue
        g=gcd(v,r)
        q,h=r//g,v//g
        mass=F(0)
        for f in fs:
            p=(-pow(f,-1,r))%r
            left=F(p,r)+a/(n*r)
            right=F(p,r)+b/(n*r)
            mass+=band_mass(n,v,left,right)
        assert mass<=F(2*r,n)+2*r/(h*(b-a))
        if q!=1 and mass:
            assert h*b>=n-q
        if g==1:
            assert mass<=2
        if q==1:
            assert mass<=len(fs)*(F(2,n)+2/(h*(b-a)))
        total+=1
    return total


def transcripts():
    refined=ROOT.joinpath('refined.txt').read_text().splitlines()
    assert len(refined)==169 and int(refined[-1].split()[0])==254
    screened=ROOT.joinpath('screen-sizes.txt').read_text().splitlines()
    assert len(screened)==135 and sum(int(l.split()[3]) for l in screened)==1555
    groups=ROOT.joinpath('ind-groups.txt').read_text().splitlines()
    assert len(groups)==1284 and max(int(l.split()[0]) for l in groups)==74
    exact=ROOT.joinpath('cover-sizes.txt').read_text().splitlines()
    counts=[list(map(int,l.split()[1:])) for l in exact if l.startswith('C ')]
    assert [r[0] for r in counts if 9<=r[0]<=42]==list(range(9,43))
    assert all(r[1]==comb(r[0]-1,4) for r in counts if 7<=r[0]<=42)
    assert sum(r[1] for r in counts if r[0]>42)==1284
    assert [l for l in exact if l.startswith('H ')]==[
        'H 7 1 2 4 5 9 12 15 18','H 8 1 3 5 7 8 10 12 14',
        'H 8 1 3 5 7 8 10 14 24','H 9 1 3 5 7 10 12 14 16']
    for line in exact:
        if line.startswith('H '):
            n,*nums=map(int,line.split()[1:])
            speeds=set(range(1,n))-set(nums[:4])|set(nums[4:])
            d=3 if n==7 else 2
            # Equality witness t=1/(dn); upper coverage is verified by the sweep.
            assert all(v%d==0 for v in speeds)
            assert min(min((v//d)%n,(-v//d)%n) for v in speeds)==1
    return refined,screened,groups,exact


def full_run(compiler):
    refined,screened,groups,exact=transcripts()
    def write(path,lines):
        with open(path,'w',encoding='ascii',newline='\n') as out:
            for l in lines:
                out.write(l+'\n')
    with tempfile.TemporaryDirectory() as directory:
        tmp=Path(directory)
        bins={}
        for name in ['reduce','refine','screen','cover']:
            exe=tmp/(name+('.exe' if os.name=='nt' else ''))
            subprocess.run([*compiler,'-O3','-std=c++17',str(ROOT/(name+'.cpp')),'-o',str(exe)],check=True)
            bins[name]=str(exe)
        def run(name,*args):
            return subprocess.check_output([bins[name],*map(str,args)],text=True).splitlines()
        first=run('reduce',11999)
        assert len(first)==2802 and first[-1].split()[0]=='3122'
        write(tmp/'first.txt',first)
        second=run('refine',tmp/'first.txt')
        assert second==refined
        write(tmp/'second.txt',second)
        print('joint-window deletion counts reproduced',flush=True)
        rows=run('screen',tmp/'second.txt',43,254)
        assert [l for l in rows if l.startswith('N ')]==screened
        rebuilt=sorted({(int(a[1]),0,*sorted(map(int,a[2:]))) for l in rows
                        if (a:=l.split())[0]=='Q'})
        assert [' '.join(map(str,r)) for r in rebuilt]==groups
        print('marked-deletion screening reproduced',flush=True)
        def exact_input():
            for n in range(7,43):
                bound=256 if n<9 else 0
                for rs in combinations(range(1,n),4):
                    yield ' '.join(map(str,(n,bound,*rs)))
            yield from groups
        write(tmp/'final.txt',exact_input())
        assert run('cover',tmp/'final.txt',7,74)==exact
        print('all exact per-size coverage rows reproduced',flush=True)


if __name__=='__main__':
    if not __debug__:
        raise RuntimeError('Run without -O so the checks are active.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full',action='store_true')
    parser.add_argument('--compiler',nargs='+',default=['c++'],help='compiler command, e.g. --compiler zig c++')
    args=parser.parse_args()
    transcripts()
    print('paired arithmetic cases',check_finite_lemmas())
    print('independent exact average samples',average_checks())
    print('first four-fold acceleration',first_four())
    small_intervals()
    if args.full:
        full_run(args.compiler)
    print('all selected checks passed')
