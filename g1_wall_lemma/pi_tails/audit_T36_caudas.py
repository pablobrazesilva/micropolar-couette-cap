"""Auditoria em Arb das ENTRADAS DE CAUDA das tabelas T3-T6 (cert_T36.py): as celulas de cauda em [-40,-20] e [30,40]
(h = 1/16) e as constantes |t| > 40 foram calculadas em double a partir das cotas rigorosas do tail_side (Lemas B, D).
Aqui cada quantidade e' recalculada em Arb e confere-se que o double usado e' >= (cotas superiores P0, P1b, rA, P0F,
rAfar) ou <= (cotas inferiores Flo, FloR) o valor Arb.  Se tudo passa, as tabelas certificadas continuam validas."""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, math, json
sys.path.insert(0, _os.path.join(_R, 'g0_constants')); sys.path.insert(0, _os.path.join(_R, 'g0_sector_lemma'))
from cl import *
from flint import arb, ctx
ctx.prec = 80
eta = '0.05'; ETA = 0.05; ETAa = arb('0.05')
TL, TR, TF = 20.0, 30.0, 40.0
bad = []; nchk = 0
def chk(name, dbl, arbv, upper=True):
    global nchk; nchk += 1
    if upper and not (dbl >= ub(arbv)): bad.append((name, dbl, ub(arbv)))
    if not upper and not (dbl <= lb(arbv)): bad.append((name, dbl, lb(arbv)))
def cell(t1, t2, side):
    # --- double, exatamente como em cert_T36.tail_vals
    T = min(abs(t1), abs(t2)); Tm = max(abs(t1), abs(t2))
    tr = tail_side(T, eta, side); tre = tail_side(T, eta, 'left' if side == 'right' else 'right')
    eB, eD = tr[3], max(tr[4], tre[4]); x = math.sqrt(T)
    P0 = 1/(2*math.pi*(2*x - 2*eD))*(1 + 1e-9); rmax = (Tm**2 + ETA**2)**0.25; rA = (rmax + eD)*(1 + 1e-9)
    P1b = rA*P0; Flo = x/(1 + eB)*(1 - 1e-9)
    # --- Arb
    Ta, Tma = arb(T), arb(Tm); xa = Ta.sqrt(); eBa, eDa = arb(eB), arb(eD)
    P0a = 1/(2*PI*(2*xa - 2*eDa)); rAa = (Tma*Tma + ETAa*ETAa).sqrt().sqrt() + eDa
    chk('P0[%g,%g]' % (t1, t2), P0, P0a); chk('rA[%g,%g]' % (t1, t2), rA, rAa); chk('P1b[%g,%g]' % (t1, t2), P1b, rAa*P0a)
    chk('Flo[%g,%g]' % (t1, t2), Flo, xa/(1 + eBa), upper=False)
for k in range(int((TF - TL)*16)): cell(-TF + k/16, -TF + (k + 1)/16, 'left')
for k in range(int((TF - TR)*16)): cell(TR + k/16, TR + (k + 1)/16, 'right')
# FloR (reflexo) nas celulas de cauda fora do compacto do G0, como em cert_T36.Flo_reflect
for k in range(int((TF - TL)*16)):
    a1, a2 = -TF + k/16, -TF + (k + 1)/16
    for (u1, u2) in ((a1, a2), (-a2, -a1)):
        T = min(abs(u1), abs(u2)); side = 'left' if -u1 <= 0 else 'right'
        if -u2 >= -TL - 1e-12 and -u1 <= TR + 1e-12: continue     # dentro do compacto: Flo certificado do G0
        tr = tail_side(max(T, 1e-9), eta, side); dbl = math.sqrt(T)/(1 + tr[3])*(1 - 1e-9)
        chk('FloR[%g,%g]' % (u1, u2), dbl, arb(T).sqrt()/(1 + arb(tr[3])), upper=False)
# constantes |t| > 40
eDF = max(tail_side(TF, eta, 'left')[4], tail_side(TF, eta, 'right')[4])
P0F = 1/(2*math.pi*(2*math.sqrt(TF) - 2*eDF)); P0Fa = 1/(2*PI*(2*arb(TF).sqrt() - 2*arb(eDF)))
chk('P0F', P0F*(1 + 1e-9), P0Fa)        # cert_T36 multiplica os termos far por (1+1e-9) no final (min com cotas uniformes)
for l in (0, 1, 4, 8, 16):
    x = TF + 16 + l + 1; d = ((x**2 + ETA**2)**0.25 + eDF)*(1 + 1e-9)
    chk('rAfar(l=%g)' % l, d, (arb(x)**2 + ETAa**2).sqrt().sqrt() + arb(eDF))
print(json.dumps(dict(checks=nchk, failures=bad[:20], nfail=len(bad)), indent=1))
