#!/usr/bin/env python3
"""Independent exact checker: read emitted TSV chains, rebuild no producer calls."""

if not __debug__:
    raise SystemExit('Do not use -O or PYTHONOPTIMIZE: this checker uses assertions.')

import argparse,collections,csv,hashlib,itertools,json,platform,sys
from fractions import Fraction
from pathlib import Path

def parity(seq):
    return 1 if sum(a>b for k,a in enumerate(seq) for b in seq[k+1:])%2==0 else -1

def exact_cover():
    pairs=list(itertools.combinations(range(7),2));position={p:k for k,p in enumerate(pairs)}
    triples=list(itertools.combinations(range(7),3));masks=[sum(1<<position[p] for p in itertools.combinations(t,2)) for t in triples]
    containing=[[k for k,m in enumerate(masks) if m>>i&1] for i in range(21)]
    answers=[]
    def search(remaining,chosen):
        if not remaining:
            answers.append(tuple(sorted(triples[i] for i in chosen)));return
        pair=(remaining&-remaining).bit_length()-1
        for triple in containing[pair]:
            mask=masks[triple]
            if mask&remaining==mask:search(remaining^mask,chosen+(triple,))
    search((1<<21)-1,())
    assert len(answers)==len(set(answers))==30
    return set(answers)

def read_tsv(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))

def determinant(A):
    a=[list(map(Fraction,r)) for r in A];det=Fraction(1)
    for i in range(len(a)):
        pivot=next(j for j in range(i,len(a)) if a[j][i])
        if pivot!=i:a[i],a[pivot]=a[pivot],a[i];det=-det
        d=a[i][i];det*=d
        for j in range(i+1,len(a)):
            factor=a[j][i]/d
            for k in range(i+1,len(a)):a[j][k]-=factor*a[i][k]
            a[j][i]=0
    return det

def run(source,out):
    data=json.loads((source/'fano_certificate.json').read_text());T=list(map(tuple,data['triples']));assert T==list(itertools.combinations(range(7),3))
    rows=read_tsv(source/'fano_planes.tsv');F=[];s=[];orientation=[]
    for n,row in enumerate(rows):
        assert int(row['plane'])==n
        F.append(tuple(tuple(map(int,t)) for t in row['seven_triples'].split()))
        s.append(int(row['label_sign']));orientation.append(int(row['orientation']))
    assert set(F)==exact_cover();assert len(F)==30
    for face in F:
        assert len(face)==7 and face==tuple(sorted(face))
        assert all(set(a)&set(b) for a,b in itertools.combinations(face,2))
        assert set().union(*map(set,face))==set(range(7))
        assert not any(all(set(hit)&set(t) for t in face) for size in (0,1,2) for hit in itertools.combinations(range(7),size))
        assert all(any(not set(a)&set(b) for a in face) for b in T if b not in face)
    base=F[0];fi={face:i for i,face in enumerate(F)}
    for p in itertools.permutations(range(7)):
        image=[tuple(sorted(p[x] for x in t)) for t in base];target=fi[tuple(sorted(image))]
        assert parity(image)==orientation[target] and parity(p)==s[target]
    chains=collections.defaultdict(dict)
    for row in read_tsv(source/'fano_cycles.tsv'):
        cycle=int(row['cycle']);coefficient=int(row['coefficient']);face=tuple(tuple(map(int,t.split(','))) for t in row['oriented_6_simplex'].split())
        assert face==tuple(sorted(face)) and len(face)==7 and len(set(face))==7
        assert all(t in T for t in face) and all(set(a)&set(b) for a,b in itertools.combinations(face,2))
        assert coefficient in (-1,1) and face not in chains[cycle]
        chains[cycle][face]=coefficient
    assert set(chains)==set(range(29))
    matrix=[]
    for n in range(29):
        c=chains[n];assert len(c)==128
        bd=collections.defaultdict(int)
        for face,coefficient in c.items():
            for i in range(7):bd[face[:i]+face[i+1:]]+=coefficient*(-1 if i%2 else 1)
        assert not any(bd.values())
        row=[c.get(face,0)*orientation[i] for i,face in enumerate(F)]
        assert sum(row[i]*s[i] for i in range(30))==0
        a,b=data['spanning_tree'][n]
        assert row==[int(i in (a,b)) for i in range(30)]
        matrix.append(row)
    reported=read_tsv(source/'fano_pairing.tsv')
    assert matrix==[[int(row[str(i)]) for i in range(30)] for row in reported]
    minor_det=determinant([r[1:] for r in matrix]);assert abs(minor_det)==1
    out.mkdir(parents=True,exist_ok=True)
    digest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.iterdir()) if p.is_file()}
    digest['certify_fano.py']=hashlib.sha256(Path(__file__).with_name('certify_fano.py').read_bytes()).hexdigest()
    digest['audit_fano.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result={'status':'PASS','independent_STS_exact_cover_count':30,'checked_ground_permutations':5040,
            'integer_cycles':29,'terms_per_cycle':128,'zero_boundaries':29,'pairing_minor_determinant':int(minor_det),
            'software':{'python':sys.version,'executable':sys.executable,'platform':platform.platform()},
            'hashes':digest,'unproved':['rational homology rank29','all-m homology rank','torsion absence','bibliographic originality']}
    (out/'independent_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.source,a.output)
