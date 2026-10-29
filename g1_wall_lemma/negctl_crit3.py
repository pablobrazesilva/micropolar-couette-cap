"""Negative control for the refined criterion (crit3arb: crit3 + Lemma F8 + route 3, no route 4), Appendix D.
At a single sigma and beta in [0, 1e-9], try all 32 parameter choices (lambda in {0,1}, theta_rho in {0.35,0.5,0.65,0.8},
variant in {(i),(ii),(ir),(iir)}) and report the best rigorous margin.  Below sigma_0' = 288 it must be <= 0.
(Written for the repository; the original control was run interactively with the same choices.)
usage: python3 negctl_crit3.py N J G sigma"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'heuristic'))
import crit3arb as C3
N, J, G, s = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
P = C3.params(float(N), float(J), float(G))
best = -9.0
for lam in (0.0, 1.0):
    for th in (0.35, 0.5, 0.65, 0.8):
        for cone in (False, True, 'ir', 'iir'):
            o = C3.CRIT(P, s, s, 0.0, 1e-9, dict(lam=lam, th_rho=th, cone=cone), 1)
            m = o.get('margin', -9.0)
            m = -9.0 if m is None or m != m else float(m)
            best = max(best, m)
print('(%s,%s,%s) sigma=%g, beta in [0,1e-9]: best rigorous margin over 32 choices = %.4f -> %s'
      % (N, J, G, s, best, 'NOT certified (as expected below sigma_0\')' if best <= 0 else 'certified'))
