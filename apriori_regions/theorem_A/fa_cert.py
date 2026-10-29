"""Certificacao Arb (python-flint) do Teorema A para a terna (1;1;0,5)  (N,J,G) = (1,1,1/2), delta_A = 0,05:
   para cada par (alpha_k, sigma_k) de fa_table.PAIRS (alpha_k < 7):  ok(alpha, sigma_k) para todo alpha in [alpha_k, 7]
   (caixas em alpha, bissecao).  Como ok e' decrescente em aR, isso da' F_A(alpha) >= sigma_k em [alpha_k, 7].
   O par (0, 3.7498) da' min_{alpha<=7} F_A >= 3.7498.  Controle: ok(0, 3.75) e' FALSO (F_A(0+) = 3.74988...).
   Tambem: enclosure de F_A(0+) (alpha=0 exato) por bissecao em aR."""
import sys, json
from flint import arb, ctx
import fa_table as FT
ctx.prec = 96
N, J, G = arb(1), arb(1), arb('0.5'); dl = arb('0.05'); pi = arb.pi()
lo_, hi_ = arb(2), arb(3)
f = lambda m: m.tan() + m.tanh()
for _ in range(80):
    m = (lo_ + hi_)/2
    if f(m).upper() < 0: lo_ = m
    elif f(m).lower() > 0: hi_ = m
    else: break
assert f(lo_).upper() < 0 and f(hi_).lower() > 0
mu0 = lo_.union(hi_)
Lam = lambda a: a**2 + pi**2*(pi**2 + a**2)/(pi**2 + 4*a**2)
Gam = lambda a: pi*mu0**2 + 2*arb(2).sqrt()*a**3
e = (G - J)/2; assert e.upper() < 0; ep = arb(0)                 # G < J: e < 0, e^+ = 0, rho = 1
def ok(a, aR):
    s = 1 - dl/Lam(a) - aR/Gam(a)
    if not s.lower() > 0: return False
    rho = 1 - 2*N*ep/(J*s)
    d = G*(1 + N)*(pi**2/4 + a**2) + 4*N
    lhs = aR*(4*N*abs(e)*J/s + 4*N**2*e**2/s**2) + dl*(J + 2*N*abs(e)/s)**2
    return (rho*J*d - lhs).lower() > 0
def notok(a, aR):                                                 # certifica a FALHA (controle negativo)
    s = 1 - dl/Lam(a) - aR/Gam(a)
    if s.upper() <= 0: return True
    d = G*(1 + N)*(pi**2/4 + a**2) + 4*N
    lhs = aR*(4*N*abs(e)*J/s + 4*N**2*e**2/s**2) + dl*(J + 2*N*abs(e)/s)**2
    return s.lower() > 0 and (J*d - lhs).upper() < 0
res = {}
for ak, sk in FT.PAIRS:
    if ak >= FT.A_CERT: continue
    st = [(ak, FT.A_CERT)]; n = 0; fails = []
    while st:
        a1, a2 = st.pop(); n += 1
        A = arb(a1).union(arb(a2))
        if ok(A, arb(str(sk))): continue
        if a2 - a1 < 1e-7: fails.append((a1, a2)); break
        m = (a1 + a2)/2; st += [(a1, m), (m, a2)]
    res['%g' % ak] = dict(sigma=sk, boxes=n, fails=fails)
    print('alpha in [%g, 7]: F_A >= %s  caixas=%d  falhas=%s' % (ak, sk, n, fails), flush=True)
# enclosure de F_A(0) por bissecao (alpha = 0 exato; F_A(0+) = F_A(0) por continuidade)
a0 = arb(0); l, h = arb('3.7'), arb('3.8')
assert ok(a0, l) and notok(a0, h)
for _ in range(40):
    m = (l + h)/2
    if ok(a0, m): l = m
    elif notok(a0, m): h = m
    else: break
neg = notok(a0, arb('3.75'))
# continuidade em alpha: ok(alpha, 3.75) falha para alpha pequeno (caixa [0, 1e-3])
neg_box = notok(arb(0).union(arb('0.001')), arb('3.75'))
print('mu0 =', mu0.str(15))
print('F_A(0) in [%s, %s];  ok(0, 3.75) falso (certificado): %s;  falso em toda alpha in [0,1e-3]: %s'
      % (l.str(12), h.str(12), neg, neg_box))
res.update(mu0=mu0.str(15), FA0=[l.str(12), h.str(12)], neg_375=neg, neg_375_box=neg_box)
json.dump(res, open('fa_cert.json', 'w'), indent=1)
