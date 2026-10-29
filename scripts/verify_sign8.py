"""Verify the 17-term certificate using only Python's standard library.
Input: data/sign8_seed.json. No discovery output or orbit/matrix data is used.
"""

if not __debug__:
    raise SystemExit('Do not use -O or PYTHONOPTIMIZE: this checker uses assertions.')

import itertools,json,time,sys,hashlib,argparse
from pathlib import Path
from collections import defaultdict
p=Path(__file__).resolve().parent
a=argparse.ArgumentParser(description=__doc__)
a.add_argument('--seed',type=Path,default=p.parent/'data/sign8_seed.json')
a.add_argument('--output',type=Path,required=True)
a=a.parse_args();t0=time.monotonic()
triples=list(itertools.combinations(range(8),3));idx={t:i for i,t in enumerate(triples)}
def sign(a):return (-1)**sum(a[i]>a[j]for i in range(len(a))for j in range(i+1,len(a)))
certificate=json.loads(a.seed.read_text());assert certificate['ground_set']==list(range(8))
seeds=[]
for r in certificate['seed_chain']:
 ts=[tuple(t)for t in r['triples']];assert len(ts)==7
 assert ts==sorted(ts) and len(set(ts))==7 and all(t in idx for t in ts)
 assert all(set(a)&set(b)for a,b in itertools.combinations(ts,2))
 assert isinstance(r['coefficient'],int)
 seeds.append((r['coefficient'],ts))
C=defaultdict(int);permutation_count=0
for g in itertools.permutations(range(8)):
 permutation_count+=1;sg=sign(g)
 for coefficient,ts in seeds:
  image=[idx[tuple(sorted(g[i]for i in t))]for t in ts]
  C[sum(1<<i for i in image)]+=coefficient*sg*sign(image)
C={m:c for m,c in C.items()if c}
D=defaultdict(int)
for m,c in C.items():
 face=[i for i in range(56)if m>>i&1]
 for j,i in enumerate(face):D[m^(1<<i)]+=(-1)**j*c
assert all(c==0 for c in D.values())
F=[]
for m,c in C.items():
 ts=[triples[i]for i in range(56)if m>>i&1];support=set().union(*map(set,ts))
 if len(support)==7 and len({a for t in ts for a in itertools.combinations(t,2)})==21:F.append((m,c))
assert len(F)==240 and {abs(c)for m,c in F}=={168}
# Check the base Fano face is a universally maximal face: inside its support
# only its seven own triples meet every line; no transversal of size <=2 exists.
base=[tuple(t)for t in certificate['seed_chain'][-1]['triples']];S=set().union(*map(set,base))
assert len(S)==7
hits=lambda A:all(set(A)&set(t)for t in base)
assert all(not hits(A)for k in range(3)for A in itertools.combinations(sorted(S),k))
assert all(not hits(t)or t in base for t in itertools.combinations(sorted(S),3))
# A triple outside S meets S in <=2 points, so cannot extend this Fano face.
# Independently reconstruct the normalized Fano sign vector on every support.
# Record label parity and simplex-orientation parity separately.
base7=list(itertools.combinations(range(7),3))
F0=tuple(t for t in base7 if (t[0]+1)^(t[1]+1)^(t[2]+1)==0)
local={}
for h in itertools.permutations(range(7)):
    im=[tuple(sorted(h[x] for x in t)) for t in F0]
    face=tuple(sorted(im));datum=(sign(h),sign(im))
    if face in local:assert local[face]==datum
    else:local[face]=datum
assert len(local)==30
normalization=[]
for omitted in range(8):
    support=[i for i in range(8) if i!=omitted]
    expected=168*(-1)**(7-omitted)
    for face,(label_sign,orientation) in local.items():
        image=[tuple(support[x] for x in t) for t in face]
        mask=sum(1<<idx[t] for t in image)
        assert C.get(mask,0)*orientation*label_sign==expected
    normalization.append(expected)
result={'status':'PASS','evidence':'EXACT_EXHAUSTIVE','seed_terms':len(seeds),'group_permutations':permutation_count,'chain_terms':len(C),'boundary_faces':len(D),'nonzero_boundary_coefficients':sum(c!=0 for c in D.values()),'fano_faces':len(F),'fano_absolute_coefficients':[168],'universal_fano_maximality_finite_conditions':True,'python':sys.version,'certificate_sha256':hashlib.sha256(a.seed.read_bytes()).hexdigest(),'normalization_by_omitted_label':normalization,'identity':'Phi(C)=-168*iota*d8(e0^...^e7)','elapsed_seconds':time.monotonic()-t0}
a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
