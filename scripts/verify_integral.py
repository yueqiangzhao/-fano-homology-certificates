#!/usr/bin/env python3
"""Exact checker for the primitive eight-label integral chain.
Regenerates the local Fano geometry without reading the supplied local-cycle data.
The only mathematical input file is data/X8_correction.tsv.
Does not compute global homology ranks or prove the all-label theorems.
"""

if not __debug__:
    raise SystemExit('Do not use -O or PYTHONOPTIMIZE: this checker uses assertions.')

from collections import Counter, defaultdict, deque
from itertools import permutations, combinations, product
from pathlib import Path
import argparse, csv, hashlib, json, re

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, default=ROOT/'verification/integral',
                    help='Directory for exact_results.json and all coordinate values')
args = parser.parse_args()
OUT = args.output.resolve()
OUT.mkdir(parents=True, exist_ok=True)

def sign_order(seq):
    return (-1) ** sum(seq[i] > seq[j] for i in range(len(seq)) for j in range(i+1,len(seq)))

def canonical_simplex(vertices):
    assert len(set(vertices)) == len(vertices)
    return tuple(sorted(vertices)), sign_order(vertices)

def boundary(chain):
    out = defaultdict(int)
    for simplex,c in chain.items():
        for i in range(len(simplex)):
            out[simplex[:i]+simplex[i+1:]] += c * (-1)**i
    return out

def is_face(simplex):
    return all(set(a)&set(b) for a,b in combinations(simplex,2))

# Transported orientation parity and label parity are different quantities.
F0 = ((0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,6),(2,3,6),(2,4,5))
planes = {}
reps = Counter()
for perm in permutations(range(7)):
    transported = tuple(tuple(sorted(perm[i] for i in triple)) for triple in F0)
    key, orientation = canonical_simplex(transported)
    values = (orientation, sign_order(perm))
    if key in planes:
        assert planes[key] == values, ('inconsistent transport', perm)
    planes[key] = values
    reps[key] += 1
assert len(planes) == 30 and set(reps.values()) == {168}
for F in planes:
    pairs=Counter(pair for t in F for pair in combinations(t,2))
    assert len(pairs)==21 and set(pairs.values()) == {1}
    assert is_face(F)
    # Complete finite cases for possible additional triples: intersections with
    # the support have <=2 labels, or are nonline triples within the support.
    for k in range(3):
        assert all(any(not(set(U)&set(line)) for line in F) for U in combinations(range(7),k))
    for B in combinations(range(7),3):
        if B not in F:
            assert any(not(set(B)&set(line)) for line in F)

# Independently regenerate local graph and 29 integral cycles, not supplied tree.
keys=sorted(planes)
edges=[(i,j) for i,j in combinations(range(30),2) if not(set(keys[i])&set(keys[j]))]
assert len(edges)==120
adj=defaultdict(list)
for i,j in edges:
    assert planes[keys[i]][1] == -planes[keys[j]][1]
    adj[i].append(j); adj[j].append(i)
assert {len(adj[i]) for i in range(30)}=={8}
seen={0}; tree=[]; q=deque([0])
while q:
    i=q.popleft()
    for j in adj[i]:
        if j not in seen:
            seen.add(j);q.append(j);tree.append((i,j))
assert len(seen)==30 and len(tree)==29
local_boundary_sizes=[]
for i,j in tree:
    F,G=keys[i],keys[j]
    antipodes=[]
    for a in F:
        miss=[b for b in G if not(set(a)&set(b))]
        assert len(miss)==1
        antipodes.append(miss[0])
    assert set(antipodes)==set(G)
    chain=Counter()
    for mask in product([0,1],repeat=7):
        simplex=tuple(antipodes[t] if mask[t] else F[t] for t in range(7))
        assert is_face(simplex)
        key,sgn=canonical_simplex(simplex)
        chain[key] += planes[F][0]*sgn*(-1)**sum(mask)
    assert len(chain)==128 and set(chain.values()) <= {-1,1}
    bd=boundary(chain)
    assert not any(bd.values())
    local_boundary_sizes.append(len(bd))
    evals={p:chain.get(p,0)*planes[p][0] for p in keys}
    assert evals[F]==evals[G]==1
    assert sum(c!=0 for c in evals.values())==2
# Root-deleted unsigned incidence minor: exact leaf elimination gives unit pivots.
remaining=set(range(30)); active=set(range(29)); pivots=[]
while active:
    leaf=next(v for v in remaining if v!=0 and sum(v in tree[e] for e in active)==1)
    edge=next(e for e in active if leaf in tree[e])
    pivots.append(1);active.remove(edge);remaining.remove(leaf)
assert remaining=={0} and pivots==[1]*29

# Read the whole literal eight-label chain, without auxiliary checker/results.
chain=Counter(); rows=0; original_coeffs=Counter(); raw_ordered=0
with (ROOT/'data/X8_correction.tsv').open(newline='') as f:
    reader=csv.DictReader(f,delimiter='\t')
    assert reader.fieldnames==['chain','coefficient','oriented_simplex_of_triples']
    for row in reader:
        assert row['chain']=='Gamma8'
        c=int(row['coefficient']); assert c
        original_coeffs[c]+=1
        triples=tuple(tuple(map(int,x.split(','))) for x in re.findall(r'\(([^()]*)\)',row['oriented_simplex_of_triples']))
        assert len(triples)==7
        assert all(len(t)==3 and len(set(t))==3 and set(t)<=set(range(8)) for t in triples)
        # Each triple is a simplicial vertex, not an internally oriented face.
        triples=tuple(tuple(sorted(t)) for t in triples)
        assert is_face(triples)
        key,sgn=canonical_simplex(triples)
        raw_ordered += (triples==key)
        chain[key]+=c*sgn
        rows+=1
assert rows==4100 and len(chain)==4100 and all(chain.values())
bd=boundary(chain)
assert len(bd)==14257 and not any(bd.values())
blocks=[];all_evals=[]
for omitted in range(8):
    T=tuple(i for i in range(8) if i!=omitted)
    vals=[]
    for F in keys:
        lifted=tuple(tuple(T[i] for i in triple) for triple in F)
        alpha=chain.get(lifted,0)*planes[F][0]
        s=planes[F][1]
        vals.append((F,alpha,s,lifted))
        all_evals.append({'omitted':omitted,'plane':[list(t) for t in lifted], 'alpha':alpha,'sign':s})
    total=sum(a*s for _,a,s,_ in vals)
    assert total==(-1)**omitted
    assert sum(a!=0 for _,a,_,_ in vals)==1
    nonzero=next({'plane':[list(t) for t in L],'alpha':a,'sign':s} for F,a,s,L in vals if a)
    blocks.append({'omitted':omitted,'signed_sum':total,'nonzero_count':1,'nonzero':nonzero})
Fstar=((1,2,3),(1,4,5),(1,6,7),(2,4,6),(2,5,7),(3,4,7),(3,5,6))
assert Fstar == tuple(tuple(i+1 for i in t) for t in F0)
assert planes[F0]==(1,1) and chain.get(Fstar,0)==0

result={
 'status':'PASS',
 'input_sha256':hashlib.sha256((ROOT/'data/X8_correction.tsv').read_bytes()).hexdigest(),
 'fano':{'permutations':5040,'planes':30,'representatives_per_plane':168,'transport_consistency':True,'all_maximality_patterns_checked':True},
 'tree':{'edges_in_graph':len(edges),'degrees':8,'tree':tree,'cycles':29,'terms_each':128,'boundary_zero':True,'boundary_faces_each':sorted(set(local_boundary_sizes)),'evaluations':'exactly u_F + u_G','root_deleted_integer_pivots':pivots},
 'Gamma8':{'raw_terms':rows,'distinct_canonical_terms':len(chain),'lex_ordered_input_rows':raw_ordered,'raw_coefficients':dict(sorted(original_coeffs.items())),'boundary_faces':len(bd),'nonzero_boundary_faces':sum(v!=0 for v in bd.values()),'fano_evaluations':len(all_evals),'blocks':blocks,'Fstar':{'sign':1,'alpha':0}},
 'scope':'Finite calculations certify the stated chain and local integral pairings only. Uniform passage to all finite sets is a separate paper proof.'
}
(OUT/'exact_results.json').write_text(json.dumps(result,indent=2)+'\n')
(OUT/'all_240_fano_evaluations.json').write_text(json.dumps(all_evals,indent=2)+'\n')
print(json.dumps(result,indent=2))
