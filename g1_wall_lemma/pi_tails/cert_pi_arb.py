"""Recertificacao INTEGRAL em Arb (python-flint) de Pi, Pi', Pi_psi, Pi'_psi na reta z(t) = t e^{-i pi/6} - eta e^{i pi/3},
eta = 0.05 (bordo de H_delta), com as CAUDAS tambem em bolas (ressalva (c) da Obs. rem:status do artigo).
  Pi   = sup|Ai(Z) Ai(wb Z)|,  Pi'  = sup max(|Ai'(Z)Ai(wb Z)|, |Ai(Z)Ai'(wb Z)|),
  Pi_psi = sup|Ai(wb Z) psi(Z)|, Pi'_psi = sup|Ai'(wb Z) psi(Z)|,   psi = Ai/F,  wb = e^{-2 pi i/3}.
Compacto: Pi, Pi' em |t| <= 40 por caixas de extremos RACIONAIS EXATOS (t_k = -40 + k/50, refinadas por bissecao exata);
          Pi_psi, Pi'_psi nas 12 800 caixas do G0 (line_0.05.json) em [-20, 30] -- contiguidade conferida.
Caudas:   identidade do Wronskiano Ai(Z)Ai(wbZ) = e^{i pi/6}/(2 pi [wb (Ai'/Ai)(wbZ) - (Ai'/Ai)(Z)]) (conferida em Arb),
          Ai'/Ai = -sqrt z + f, |f| <= eps_D (Lema D) em Z e em wbZ = conj Z(-t);  wb sqrt(wb Z) = -sqrt Z  =>
          |denominador| >= 2|Z|^{1/2} - 2 eps_D;  |Ai'/Ai| <= |Z|^{1/2} + eps_D;  |F| >= |Z|^{1/2}/(1+eps_B) (Lema B).
          As cotas sao DECRESCENTES em r = |Z|^{1/2} (r >= sqrt|t| >= sqrt T), logo avaliadas em r = sqrt T, tudo em Arb.
uso: python3 cert_pi_arb.py   -> cert_pi_arb.json"""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, json, math, time
from fractions import Fraction as Fr
sys.path.insert(0, _os.path.join(_R, 'g0_constants')); sys.path.insert(0, _os.path.join(_R, 'g0_sector_lemma'))
from cl import *                          # tail_side (Arb), eps_D, ub, lb, PI
from flint import arb, acb, ctx
ctx.prec = 80
ETA = arb('0.05'); eta = '0.05'
EM = acb(arb(3).sqrt()/2, -arb(1)/2); E3 = acb(arb(1)/2, arb(3).sqrt()/2); WB = acb(-arb(1)/2, -arb(3).sqrt()/2)
W6 = acb(arb(3).sqrt()/2, arb(1)/2)                                   # e^{i pi/6}
T0 = time.time()
def A(q): return arb(q.numerator)/q.denominator                     # exact rational -> ball containing it
def zball(t1, t2):
    """acb (retangulo) contendo z(t) para todo t no intervalo EXATO [t1,t2] (racionais)"""
    tc = (A(t1) + A(t2))/2; r = (A(t2) - A(t1))/2
    zc = acb(tc)*EM - ETA*E3
    m = acb(arb(zc.real.mid()), arb(zc.imag.mid()))
    rad = ub(r) + ub(abs(zc - m))                                    # |z(t) - m| <= |t - tc| + |zc - m|
    return acb(arb(zc.real.mid(), rad), arb(zc.imag.mid(), rad))
# ---- (0) a identidade do Wronskiano, conferida em Arb em 3 pontos
for zz in (acb('0.3', '0.7'), acb('-2.1', '1.3'), acb('5', '-0.4')):
    a, ap, _, _ = zz.airy(); b, bp, _, _ = (WB*zz).airy()
    lhs = a*b; rhs = W6/(2*PI*(WB*bp/b - ap/a))
    assert abs(lhs - rhs).upper() < 1e-20, (zz, lhs, rhs)
# ---- (1) Pi, Pi', Pi'b no compacto |t| <= 40
T = 40; n = 4000; h = Fr(1, 50)
_old = json.load(open(_os.path.join(_CERT, 'g1_wall_lemma', 'cert_pi_0.05.json'))); U0, U1 = _old['Pi'], _old['Pi1']
P0 = P1 = P1b = arb(0); at0 = None; nbox = 0
for k in range(n):
    st = [(Fr(-T) + k*h, Fr(-T) + (k + 1)*h)]
    while st:
        a, b = st.pop()
        Zb = zball(a, b); ai, aip, _, _ = Zb.airy(); bi, bip, _, _ = (WB*Zb).airy(); nbox += 1
        p0, p1, p1b = ub(abs(ai*bi)), ub(abs(aip*bi)), ub(abs(ai*bip))
        if b - a > h/4096 and (p0 > U0 - 2e-6 or p1 > U1 - 2e-6 or p1b > U1 - 2e-6):
            m = (a + b)/2; st += [(a, m), (m, b)]; continue
        if p0 > ub(P0): P0, at0 = arb(p0), (float(a), float(b))
        P1 = arb(max(ub(P1), p1)); P1b = arb(max(ub(P1b), p1b))
t_comp = time.time() - T0
# ---- (2) caudas |t| >= 40, em Arb
def side_eps(Tv, side):
    r = tail_side(Tv, eta, side)          # (ulo, utlo, Fup, eB, eD) -- rigorosos (Arb), devolvidos como floats outward
    return arb(r[3]), arb(r[4])
eBr, eDr = side_eps(T, 'right'); eBl, eDl = side_eps(T, 'left')
eD = eDr if ub(eDr) >= ub(eDl) else eDl
x = arb(T).sqrt()
assert lb(x - eD) > 0
PiT = 1/(2*PI*(2*x - 2*eD)); Pi1T = (x + eD)*PiT
Pi = max(ub(P0), ub(PiT)); Pi1 = max(ub(P1), ub(P1b), ub(Pi1T))
# ---- (3) Pi_psi, Pi'_psi: compacto nas caixas do G0 + caudas
d = json.load(open(_os.path.join(_CERT, 'g0_constants', 'line_0.05.json'))); boxes = d['boxes']
ts = sorted((bx[0], bx[1]) for bx in boxes)
assert ts[0][0] == -20.0 and ts[-1][1] == 30.0 and all(ts[i][1] == ts[i + 1][0] for i in range(len(ts) - 1)), 'caixas G0 nao contiguas'
Q0 = Q1 = 0.0; atq0 = atq1 = None
_oq = json.load(open(_os.path.join(_CERT, 'g1_wall_lemma', 'cert_pipsi.json'))); V0, V1 = _oq['Pi_psi'], _oq['Pi1_psi']
for (t1, t2, Flo, Fup, *_r) in boxes:
    st = [(Fr(t1), Fr(t2))]
    while st:                                  # sub-bolas da caixa do G0; Flo vale em toda a caixa, logo nas sub-bolas
        a, b = st.pop()
        Zb = zball(a, b); ai, aip, _, _ = Zb.airy(); bi, bip, _, _ = (WB*Zb).airy()
        q0 = ub(abs(ai*bi)/arb(Flo)); q1 = ub(abs(ai*bip)/arb(Flo))
        if (b - a) > Fr(1, 256*4096) and (q0 > V0 - 2e-6 or q1 > V1 - 2e-6):
            m = (a + b)/2; st += [(a, m), (m, b)]; continue
        if q0 > Q0: Q0, atq0 = q0, (float(a), float(b))
        if q1 > Q1: Q1, atq1 = q1, (float(a), float(b))
def tail_psi(Tv, here, refl):
    eB, eDh = side_eps(Tv, here); _, eDr_ = side_eps(Tv, refl)
    e = eDh if ub(eDh) >= ub(eDr_) else eDr_
    xx = arb(Tv).sqrt(); assert lb(xx - e) > 0 and ub(eB) < 1
    PT = 1/(2*PI*(2*xx - 2*e))
    return ub(PT*(1 + eB)/xx), ub((xx + e)*PT*(1 + eB)/xx)
tl = tail_psi(20, 'left', 'right'); trr = tail_psi(30, 'right', 'left')
Pipsi = max(Q0, tl[0], trr[0]); Pi1psi = max(Q1, tl[1], trr[1])
# ---- (4) comparacao com os valores usados no G1 (cert_pi_0.05.json, cert_pipsi.json)
old_pi = json.load(open(_os.path.join(_CERT, 'g1_wall_lemma', 'cert_pi_0.05.json'))); old_ps = json.load(open(_os.path.join(_CERT, 'g1_wall_lemma', 'cert_pipsi.json')))
used = dict(Pi=old_pi['Pi'], Pi1=old_pi['Pi1'], Pi_psi=old_ps['Pi_psi'], Pi1_psi=old_ps['Pi1_psi'])
new = dict(Pi=Pi, Pi1=Pi1, Pi_psi=Pipsi, Pi1_psi=Pi1psi)
out = dict(eta=eta, T_Pi=T, h='1/50 (exato), bissecao exata', boxes_Pi=nbox, Pi_box=ub(P0), Pi_box_at=at0, Pi1_box=ub(P1), Pi1b_box=ub(P1b),
           eps_D_40=ub(eD), tail_Pi=ub(PiT), tail_Pi1=ub(Pi1T), Pi_psi_box=Q0, Pi_psi_at=atq0, Pi1_psi_box=Q1, Pi1_psi_at=atq1,
           tail_psi_left_T20=tl, tail_psi_right_T30=trr, new=new, used_in_G1=used,
           all_used_are_upper_bounds=all(used[k] >= new[k] for k in new), seconds=time.time() - T0, seconds_compact_Pi=t_comp)
print(json.dumps(out, indent=1)); json.dump(out, open(_wout('g1_wall_lemma', 'cert_pi_arb.json'), 'w'), indent=1)
