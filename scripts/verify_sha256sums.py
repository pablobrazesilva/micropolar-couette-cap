#!/usr/bin/env python3
"""Check SHA256SUMS (same format as `sha256sum -c`), Python standard library only."""
import hashlib, os, sys
R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
bad = n = 0
for line in open(os.path.join(R, 'SHA256SUMS')):
    if not line.strip(): continue
    h, f = line.rstrip('\n').split(None, 1); f = f.lstrip('*'); n += 1
    p = os.path.join(R, f)
    if not os.path.exists(p) or hashlib.sha256(open(p, 'rb').read()).hexdigest() != h:
        bad += 1; print('  MISMATCH or MISSING:', f)
print('  SHA256SUMS: %d files, %d mismatches' % (n, bad))
sys.exit(1 if bad else 0)
