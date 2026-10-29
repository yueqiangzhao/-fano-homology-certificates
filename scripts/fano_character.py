#!/usr/bin/env python3
"""All 15 class traces of the actual 30-plane permutation action.
The standard two-row character identity chi_(4,3)=M_3-M_2 is compared to
fixed-plane counts, giving a representation-theoretic bridge separately
from the old Lefschetz trace calculation.
"""

if not __debug__:
    raise SystemExit('Do not use -O or PYTHONOPTIMIZE: this checker uses assertions.')

import argparse,collections,csv,itertools,json,math
from pathlib import Path
from fractions import Fraction

def parts(n,bound=None):
    if not n:yield ();return
    for k in range(min(n,bound or n),0,-1):
        for tail in parts(n-k,k):yield (k,)+tail

def run(source,output):
    d=json.loads(source.read_text());T=d['triples'];planes=[tuple(tuple(T[i]) for i in F['face']) for F in d['planes']];rows=[]
    for mu in parts(7):
        perm=list(range(7));offset=0
        for k in mu:
            for j in range(k):perm[offset+j]=offset+(j+1)%k
            offset+=k
        fixed=sum(tuple(sorted(tuple(sorted(perm[x] for x in t)) for t in F))==F for F in planes)
        subset_traces=[]
        for k in (2,3):
            subset_traces.append(sum(tuple(sorted(perm[x] for x in A))==A for A in itertools.combinations(range(7),k)))
        sg=(-1)**(7-len(mu));chi=subset_traces[1]-subset_traces[0]
        assert fixed==1+sg+chi+sg*chi
        z=math.prod(k**a*math.factorial(a) for k,a in collections.Counter(mu).items())
        rows.append([mu,z,fixed,sg,chi,sg*chi,fixed-sg])
    for a,b in itertools.product(range(2,7),repeat=2):
        val=sum(Fraction(r[a]*r[b],r[1]) for r in rows)
        if a==b and a in (3,4,5):assert val==1
        if a!=b and a in (3,4,5) and b in (3,4,5):assert val==0
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['cycle_type','z','fixed_Fano_planes','sign','chi_4_3','chi_2_2_2_1','H6_character'])
        for row in rows:w.writerow(['+'.join(map(str,row[0]))]+row[1:])
    print('PASS: all 15 actual Fano traces = 1+sign+chi_4_3+chi_2_2_2_1; three nontrivial constituents have exact orthonormal characters')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.source,a.output)
