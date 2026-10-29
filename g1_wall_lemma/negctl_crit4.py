"""Negative control for the criterion with route 4 (crit4arb: crit3 + Lemma F8 + routes 3 and 4), Appendix D.
At a single sigma and beta in [0, 1e-9], try all 32 parameter choices (lambda in {0,1}, theta_rho in {0.35,0.5,0.65,0.8},
variant in {(i),(ii),(ir),(iir)}) and report the best rigorous margin.  Below the new sigma_0' it must be <= 0 at small beta.
(Same 32 choices as negctl_crit3.py.)
usage: python3 negctl_crit4.py N J G sigma [beta1 beta2]   (default beta in [0, 1e-9]); the box is [sigma, 1.005 sigma] x [beta1, beta2]"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'heuristic'))
import crit4arb as C3
N, J, G, s = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
b1, b2 = (float(sys.argv[5]), float(sys.argv[6])) if len(sys.argv) > 6 else (0.0, 1e-9)
whys = set()
P = C3.params(float(N), float(J), float(G))
best = -9.0
for lam in (0.0, 1.0):
    for th in (0.35, 0.5, 0.65, 0.8):
        for cone in (False, True, 'ir', 'iir'):
            o = C3.CRIT(P, s, s*1.005, b1, b2, dict(lam=lam, th_rho=th, cone=cone), 1); whys.add(o.get('why'))
            m = o.get('margin', -9.0)
            m = -9.0 if m is None or m != m else float(m)
            best = max(best, m)
print('(%s,%s,%s) sigma in [%g,%g], beta in [%g,%g]: best rigorous margin over 32 choices = %.4f -> %s  (reasons: %s)'
      % (N, J, G, s, s*1.005, b1, b2, best, 'NOT certified (as expected below sigma_0\')' if best <= 0 else 'certified', ','.join(sorted(map(str, whys)))))
