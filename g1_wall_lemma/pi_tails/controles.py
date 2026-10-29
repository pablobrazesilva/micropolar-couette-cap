"""Controles para cert_pi_arb.py (Arb): (1) valores pontuais (cotas INFERIORES do sup) perto dos maximos;
(2) caudas: valores pontuais em |t| >= T ficam abaixo da cota de cauda (consistencia), em t = +-40, +-60, +-200, +-2000;
(3) negativo: a cota de cauda NAO serve com T pequeno (T = 2: eps_D >= sqrt T ou cota > Pi);
(4) negativo: um alvo abaixo do valor pontual (Pi_alvo = valor pontual - 1e-7) nao e' certificavel (o refinamento
    sempre encontra caixa com cota > alvo)."""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, json
from fractions import Fraction as Fr
sys.path.insert(0, _os.path.join(_R, 'g0_constants')); sys.path.insert(0, _os.path.join(_R, 'g0_sector_lemma'))
from cl import *
from flint import arb, acb, ctx
ctx.prec = 80
ETA = arb('0.05'); EM = acb(arb(3).sqrt()/2, -arb(1)/2); E3 = acb(arb(1)/2, arb(3).sqrt()/2); WB = acb(-arb(1)/2, -arb(3).sqrt()/2)
def z(t): return acb(arb(t))*EM - ETA*E3
def vals(t):
    Z = z(t); a, ap, _, _ = Z.airy(); b, bp, _, _ = (WB*Z).airy(); F, _, _ = point(str(t), '0.05', 80)
    return abs(a*b), abs(ap*b), abs(a*bp), abs(b*a/F), abs(bp*a/F)
d = json.load(open(_os.path.join(_CERT, 'g1_wall_lemma', 'cert_pi_arb.json')))
# (1) cotas inferiores pontuais
import numpy as np
best = [0]*5; arg = [None]*5
for t in list(np.arange(-0.02, 0.02, 0.0005)) + list(np.arange(0.70, 0.78, 0.001)) + list(np.arange(-0.75, -0.65, 0.001)) + list(np.arange(-1.40, -1.30, 0.001)):
    v = vals(repr(float(t)))
    for i in range(5):
        if lb(v[i]) > best[i]: best[i], arg[i] = lb(v[i]), float(t)
names = ['Pi', 'Pi1(Ai\')', 'Pi1b(Ai\'wb)', 'Pi_psi', 'Pi1_psi']
cert = [d['new']['Pi'], d['Pi1_box'], d['Pi1b_box'], d['new']['Pi_psi'], d['new']['Pi1_psi']]
for i in range(5): print('(1) %-12s cota inferior pontual %.10f (t=%.4f)   cota superior certificada %.10f   folga %.2e' % (names[i], best[i], arg[i], cert[i], cert[i] - best[i]))
# (2) caudas
for t in (40, -40, 60, -60, 200, -200, 2000, -2000):
    v = vals(str(t))
    print('(2) t=%6d |AiAi|=%.5f <= %.5f   |Ai\'Ai|=%.5f,|AiAi\'|=%.5f <= %.5f' % (t, ub(v[0]), d['tail_Pi'], ub(v[1]), ub(v[2]), d['tail_Pi1']),
          'OK' if ub(v[0]) <= d['tail_Pi'] and max(ub(v[1]), ub(v[2])) <= d['tail_Pi1'] else 'FALHA')
for t, side in ((-20, 'L'), (-25, 'L'), (30, 'R'), (40, 'R'), (-200, 'L'), (200, 'R')):
    v = vals(str(t)); bnd = d['tail_psi_left_T20'] if side == 'L' else d['tail_psi_right_T30']
    print('(2) t=%5d |Ai(wb)psi|=%.5f <= %.5f  |Ai\'(wb)psi|=%.5f <= %.5f' % (t, ub(v[3]), bnd[0], ub(v[4]), bnd[1]), 'OK' if ub(v[3]) <= bnd[0] and ub(v[4]) <= bnd[1] else 'FALHA')
# (3) T pequeno
for Tv in (2, 5, 10):
    try:
        r = tail_side(Tv, '0.05', 'right'); l = tail_side(Tv, '0.05', 'left')
    except AssertionError:
        print('(3) T=%d: hipoteses dos Lemas B/D falham (eps_B nao < 1) -> cota de cauda inaplicavel (esperado)' % Tv); continue
    eD = max(r[4], l[4]); x = arb(Tv).sqrt()
    if lb(x) <= eD: print('(3) T=%d: eps_D=%.3f >= sqrt T -> cota de cauda inaplicavel (esperado)' % (Tv, eD)); continue
    PiT = 1/(2*PI*(2*x - 2*eD)); print('(3) T=%d: cauda Pi=%.4f, Pi1=%.4f  (vs Pi=%.4f, Pi1=%.4f)' % (Tv, ub(PiT), ub((x + eD)*PiT), d['new']['Pi'], d['new']['Pi1']))
# (4) alvo abaixo do valor pontual: refinamento ate h/2^14 perto do argmax nao certifica
alvo = best[0] - 1e-7; a, b = Fr(-1, 400), Fr(1, 400); ok = True; st = [(a, b)]; n = 0
def A(q): return arb(q.numerator)/q.denominator
while st and n < 4000:
    a, b = st.pop(); n += 1
    tc = (A(a) + A(b))/2; rr = (A(b) - A(a))/2; zc = acb(tc)*EM - ETA*E3; m = acb(arb(zc.real.mid()), arb(zc.imag.mid())); rad = ub(rr) + ub(abs(zc - m))
    Zb = acb(arb(zc.real.mid(), rad), arb(zc.imag.mid(), rad)); ai, _, _, _ = Zb.airy(); bi, _, _, _ = (WB*Zb).airy()
    if ub(abs(ai*bi)) > alvo:
        if b - a < Fr(1, 50*2**14): ok = False; break
        mm = (a + b)/2; st += [(a, mm), (mm, b)]
print('(4) alvo Pi = %.10f (abaixo do valor pontual): %s' % (alvo, 'NAO certificavel (esperado)' if not ok else 'certificado?! (FALHA do controle)'))
