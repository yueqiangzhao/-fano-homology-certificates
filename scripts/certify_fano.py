#!/usr/bin/env python3
"""Exact finite certificate for a 29-element integral cycle family in X_7.

Program objects: a simplicial vertex is a sorted ground-label triple, an
oriented 6-simplex is a sorted 7-tuple of such vertices. Chain coefficients
are integers. Every emitted cycle is the join of seven vertex differences.
Only finite combinatorics is certified here; rational completeness uses the
separately stated rank dim H_6(X_7;Q)=29.
"""

if not __debug__:
    raise SystemExit('Do not use -O or PYTHONOPTIMIZE: this checker uses assertions.')

import argparse,collections,csv,itertools,json,pathlib

def sign(xs):
    return (-1)**sum(xs[i]>xs[j] for i in range(len(xs)) for j in range(i+1,len(xs)))

def matchings(xs):
    if not xs:
        yield ()
        return
    a=xs[0]
    for k in range(1,len(xs)):
        b=xs[k]
        for rest in matchings(xs[1:k]+xs[k+1:]):
            yield ((a,b),)+rest

def run(out):
    out.mkdir(parents=True,exist_ok=True)
    labels=tuple(range(7));triples=tuple(itertools.combinations(labels,3));ti={t:i for i,t in enumerate(triples)}
    base=tuple(i for i,t in enumerate(triples) if (t[0]+1)^(t[1]+1)^(t[2]+1)==0)
    planes={}
    for p in itertools.permutations(labels):
        image=tuple(ti[tuple(sorted(p[x] for x in triples[a]))] for a in base)
        F=tuple(sorted(image));datum={'transport':p,'orientation':sign(image),'label_sign':sign(p)}
        if F in planes:
            assert all(planes[F][k]==datum[k] for k in ('orientation','label_sign'))
        else:planes[F]=datum
    systematic=set()
    for pairs in matchings(tuple(range(1,7))):
        common=[tuple(sorted((0,)+pair)) for pair in pairs]
        for parity in (0,1):
            lines=common+[tuple(sorted(pairs[k][bits[k]] for k in range(3)))
                          for bits in itertools.product((0,1),repeat=3) if sum(bits)%2==parity]
            systematic.add(tuple(sorted(ti[t] for t in lines)))
    assert len(planes)==30 and set(planes)==systematic
    faces=tuple(sorted(planes));index={F:i for i,F in enumerate(faces)}
    for F in faces:
        count=collections.Counter(pair for i in F for pair in itertools.combinations(triples[i],2))
        assert count==collections.Counter(itertools.combinations(labels,2))
        assert all(set(triples[a])&set(triples[b]) for a,b in itertools.combinations(F,2))
        assert not any(all(set(hit)&set(triples[a]) for a in F) for n in (0,1,2) for hit in itertools.combinations(labels,n))
        assert all(any(not set(t)&set(triples[a]) for a in F) for t in triples if ti[t] not in F)
    adj=[[] for _ in faces];edges={}
    for i,j in itertools.combinations(range(30),2):
        F,G=faces[i],faces[j]
        assert len(set(F)&set(G)) in (0,1,3)
        if set(F)&set(G):continue
        pairs=[]
        for a in F:
            disjoint=[b for b in G if not set(triples[a])&set(triples[b])]
            assert len(disjoint)==1;pairs.append((a,disjoint[0]))
        assert len({b for a,b in pairs})==7
        assert all(set(triples[a])&set(triples[b]) for k,P in enumerate(pairs) for Q in pairs[k+1:] for a in P for b in Q)
        adj[i].append(j);adj[j].append(i);edges[i,j]=pairs
    assert len(edges)==120 and all(len(ns)==8 for ns in adj)
    seen={0};queue=collections.deque([0]);tree=[]
    while queue:
        a=queue.popleft()
        for b in sorted(adj[a]):
            if b not in seen:
                seen.add(b);queue.append(b);tree.append((a,b))
    assert len(seen)==30 and len(tree)==29
    chain_rows=[];matrix=[];cycle_data=[]
    for cycle_index,(a,b) in enumerate(tree):
        i,j=sorted((a,b));pairs=edges[i,j];c={}
        for bits in itertools.product((0,1),repeat=7):
            vertices=[pairs[k][bits[k]] for k in range(7)]
            face=tuple(sorted(vertices));assert face not in c
            c[face]=(-1)**sum(bits)*sign(vertices)
        bd=collections.Counter()
        for face,coefficient in c.items():
            for k in range(7):bd[face[:k]+face[k+1:]]+=coefficient*(-1)**k
        assert len(c)==128 and not any(bd.values())
        row=[c.get(F,0)*planes[F]['orientation'] for F in faces]
        assert row==[int(k==a or k==b) for k in range(30)]
        assert planes[faces[a]]['label_sign']==-planes[faces[b]]['label_sign']
        matrix.append(row)
        cycle_data.append({'edge':[a,b],'pairs':[list(P) for P in pairs]})
        for face,coefficient in sorted(c.items()):
            chain_rows.append([cycle_index,coefficient,' '.join(','.join(map(str,triples[v])) for v in face)])
    # Tree rows: deletion of root column gives an integer matrix with determinant +/-1.
    # Verify by stripping leaf columns with exact integer pivots; no numeric rank routine.
    remaining_vertices=set(range(30));remaining_edges=set(range(29));pivots=[]
    while remaining_edges:
        leaf=next(v for v in sorted(remaining_vertices-{0}) if sum(v in tree[e] for e in remaining_edges)==1)
        edge=next(e for e in remaining_edges if leaf in tree[e]);assert matrix[edge][leaf]==1
        pivots.append([edge,leaf]);remaining_edges.remove(edge);remaining_vertices.remove(leaf)
    assert remaining_vertices=={0}
    with (out/'fano_planes.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['plane','orientation','label_sign','seven_triples','transport'])
        for i,F in enumerate(faces):
            d=planes[F];w.writerow([i,d['orientation'],d['label_sign'],' '.join(''.join(map(str,triples[v])) for v in F),','.join(map(str,d['transport']))])
    with (out/'fano_cycles.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['cycle','coefficient','oriented_6_simplex']);w.writerows(chain_rows)
    with (out/'fano_pairing.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['cycle']+list(range(30)))
        for i,row in enumerate(matrix):w.writerow([i]+row)
    cert={'ground_labels':labels,'triples':triples,'base_plane':base,
          'planes':[{'face':F,**planes[F]} for F in faces],
          'disjointness_edges':[list(e) for e in edges],
          'spanning_tree':tree,'cycles':cycle_data,'pairing':matrix,
          'unimodular_leaf_pivots':pivots,
          'classification':'EXACT_EXHAUSTIVE_FINITE_CERTIFICATE',
          'claims':['30 planes exhaust all STS(7) by the point-0 matching construction',
                    'all 5040 transports have consistent orientations and label parity',
                    'each plane is maximal in X7 and has no <=2 transversal',
                    '120 disjoint-plane pairs give seven disjointness axes',
                    '29 explicit integral 128-term cycles are closed',
                    'their Fano coefficient pairing has a unimodular 29x29 minor'],
          'not_proved_by_program':['rational homology rank29','all-m rank formula','general homology theorem']}
    (out/'fano_certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    summary={'planes':30,'disjoint_pairs':120,'tree_edges':29,'terms_per_cycle':128,
             'integral_cycle_terms':len(chain_rows),'integer_pairing_rank':29,
             'minor_absolute_determinant':1,'status':'PASS'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=pathlib.Path,default=pathlib.Path('work/fano/verified'));a=p.parse_args();run(a.output)
