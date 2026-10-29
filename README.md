# Reproducible finite certificates for the manuscript

**Title:** Fano cycles and integral homology lattices for intersecting triples
**Author:** Yueqiang Zhao
**Journal target:** Homology, Homotopy and Applications

This supplement contains literal data and scripts for the exact finite
calculations referenced in the manuscript. All original `data/` and `scripts/`
files are byte-for-byte copied from the fully checked October 29 revision.
The summary of the prior replay is included in `verification/`.

## Reproduce

Run from the root of this extracted supplement:

```bash
python3 scripts/replay.py --output verification/new_run
```

Python's standard library is the only requirement. **Do not use `python -O`**
or set `PYTHONOPTIMIZE`, because verification uses assertions.

The replay regenerates the 30 local Fano planes, 29 local tree cycles,
checks the local character values, reconstructs the alternating eight-label
cycle from a 17-term seed and verifies its boundary and coefficients, and
checks the 4,100-term primitive integral cycle (including 14,257 boundary
faces and all 240 Fano coordinate evaluations).

The supplied report states PASS for the finite certificates previously
replayed on 29 October 2026; this is not a claim of formal verification of
the global mathematical proofs.

## Key data files

- `data/fano/` -- local Fano facets, cycles, coefficient pairings,
  and character table.
- `data/sign8_seed.json` -- alternating eight-label seed.
- `data/X8_correction.tsv` -- primitive 4,100-term integer chain.

## Limits of this certificate

The scripts check specific finite chains and arithmetic identities.
They do not reprove the external rank theorem, establish the uniform
functorial and lattice claims without the printed arguments, determine
torsion in integral homology, or compute the full homotopy type.

`SHA256SUMS` lists the hashes for all immutable data and program files.
Verify with `sha256sum -c SHA256SUMS` (from this directory).