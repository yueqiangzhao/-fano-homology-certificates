#!/usr/bin/env python3
"""Replay the rational and integral finite certificates using only the standard library.

This verifies finite chains and coefficient data, not the external rank theorem,
the all-label arguments, torsion-freeness, or the full homotopy type.
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def main() -> None:
    if not __debug__:
        raise SystemExit("Do not use -O or PYTHONOPTIMIZE: these checkers use assertions.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'verification/replay')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    jobs = [
        ('fano_generate', 'certify_fano.py', ['--output', str(out/'fano')]),
        ('fano_separate_checker', 'audit_fano.py',
         ['--source', str(out/'fano'), '--output', str(out/'independent')]),
        ('fano_character', 'fano_character.py',
         ['--source', str(out/'fano/fano_certificate.json'),
          '--output', str(out/'fano/fano_character.tsv')]),
        ('sign8', 'verify_sign8.py',
         ['--seed', str(ROOT/'data/sign8_seed.json'), '--output', str(out/'sign8.json')]),
        ('integral', 'verify_integral.py', ['--output', str(out/'integral')]),
    ]
    statuses = {}
    for name, script, options in jobs:
        command = [sys.executable, str(ROOT/'scripts'/script), *options]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        (out/(name+'.log')).write_text(result.stdout+result.stderr)
        statuses[name] = 'PASS' if result.returncode == 0 else 'FAIL'
        print(name, statuses[name], flush=True)
        if result.returncode:
            print(result.stdout+result.stderr)
            raise SystemExit(result.returncode)
    for path in sorted((ROOT/'data/fano').iterdir()):
        if path.is_file() and path.read_bytes() != (out/'fano'/path.name).read_bytes():
            raise SystemExit('Regenerated data differ from supplied data: '+path.name)
    report = {
        'status': 'PASS', 'checks': statuses,
        'regenerated_fano_files_match': True,
        'python': sys.version,
        'scope': ('Exact local Fano chains, integral unit pairings, local character, '
                  'alternating eight-label cycle, and primitive integral eight-label chain. '
                  'The external rank theorem and all-label proofs are not computer-verified.')
    }
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
