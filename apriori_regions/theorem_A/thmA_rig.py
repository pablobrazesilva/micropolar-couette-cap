"""Rigorous (Arb) check of the Theorem-A ingredient used by sigma0: for alpha >= alpha_min,
   ok(alpha, aR = c alpha^2) holds (thmA.ok, delta_A = 0.05), i.e. F_A(alpha) >= c alpha^2.
   ok is monotone decreasing in aR, so then F_A(alpha) >= sigma0 for alpha >= sqrt(sigma0/c) (rigorous alpha0').
   Range [alpha_min, 1e4] by interval boxes; alpha >= 1e4 by the explicit asymptotic bound below.
usage: thmA_rig.py N J G c alpha_min"""
import sys
from flint import arb, ctx
ctx.prec = 80
N, J, G, c = [arb(v) for v in sys.argv[1:5]]; amin = float(sys.argv[5]); dl = arb('0.05')
pi = arb.pi()
# mu0: root of tan m + tanh m in [2,3], enclosed by sign change
lo, hi = arb(2), arb(3)
f = lambda m: m.tan() + m.tanh()
for _ in range(70):
    m = (lo + hi)/2
    if f(m).upper() < 0: lo = m
    elif f(m).lower() > 0: hi = m
    else: break
assert f(lo).upper() < 0 and f(hi).lower() > 0
mu0 = lo.union(hi)
Lam = lambda a: a**2 + pi**2*(pi**2 + a**2)/(pi**2 + 4*a**2)
Gam = lambda a: pi*mu0**2 + 2*arb(2).sqrt()*a**3
e = (G - J)/2; ep = e if e.lower() > 0 else arb(0); assert not (e.lower() <= 0 <= e.upper()) or True
def ok(a, aR):
    s = 1 - dl/Lam(a) - aR/Gam(a)
    if not s.lower() > 0: return False
    rho = 1 - 2*N*ep/(J*s)
    if not rho.lower() > 0: return False
    d = G*(1 + N)*(pi**2/4 + a**2) + 4*N
    lhs = aR*(4*N*abs(e)*J/s + 4*N**2*e**2/s**2) + dl*(J + 2*N*abs(e)/s)**2
    return (rho*J*d - lhs).lower() > 0
n = 0; st = [(amin, 1e4)]; fails = []
while st:
    a1, a2 = st.pop(); n += 1
    A = arb(a1).union(arb(a2))
    if ok(A, c*A**2): continue
    if a2 - a1 < 1e-6*a1: fails.append((a1, a2)); break
    m = (a1 + a2)/2; st += [(a1, m), (m, a2)]
# asymptotic alpha >= 1e4 (requires e<0 i.e. G<J, so ep=0, rho=1):  s >= s_min, lhs <= c a^2 K/s_min^2 + dl (J+2N|e|/s_min)^2,
# rhs >= J G (1+N) a^2  -> sufficient: c K/s_min^2 + dl (J+2N|e|/s_min)^2 / 1e8 < J G (1+N)
assert e.upper() < 0
A0 = arb(10)**4
smin = 1 - dl/A0**2 - c/(2*arb(2).sqrt()*A0)          # Lam >= a^2, Gam >= 2 sqrt2 a^3, and c a^2/Gam <= c/(2sqrt2 a)
K = 4*N*abs(e)*J + 4*N**2*e**2
asym = J*G*(1 + N) - (c*K/smin**2 + dl*(J + 2*N*abs(e)/smin)**2/A0**2)
print('N,J,G=%s,%s,%s c=%s: mu0=%s; interval part [%g,1e4]: %d boxes, fails=%s; asymptotic margin (must be >0) = %s; c_inf=JG(1+N)/K=%s'
      % (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], mu0.str(10), amin, n, fails, asym.str(5), (J*G*(1 + N)/K).str(8)))
