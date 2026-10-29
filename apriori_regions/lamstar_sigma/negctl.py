"""Controles negativos.
 (a) para cada faixa (s_lo, s_hi, L*) da tabela: com L*' = L* - 1 o criterio (avaliado em PONTOS, bolas de raio ~1e-29)
     falha em algum ponto da nova fronteira {Re Lambda = L*', Im Lambda in [-L*', S/2]} U {Im Lambda = -L*', Re in [-delta', L*']},
     com sigma = s_lo e beta = 0.  Registra o pior rho (>= 1 = falha) e o ponto.
 (b) cauda: ann_eval com P0 = 30 falha em sigma = 288 (com P0 = 100 passa).
 uso: python3 negctl.py tabela.txt"""
import sys, json, math
import ap_arb as AA
from flint import arb
rows = []
for line in open(sys.argv[1]):
    line = line.split('#')[0].split()
    if len(line) >= 3: rows.append((float(line[0]), float(line[1]), float(line[2])))
out = []
import fa_table as FT
from scan_pts import pts
for s_lo, s_hi, L in rows:
    Lp = L - 1; S = 2*s_lo**(1/3); worst = (0, None, None)
    bmax = FT.alpha_max(s_lo*1.0001)*s_lo**(-1/3)
    for b in [0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.6, 1, bmax]:
        if b > bmax: continue
        for Lm in pts(Lp, S):
            x, y = Lm.real, Lm.imag
            if x < -0.0396: x = -0.0396                       # dentro da regiao Re Lambda >= -delta'
            r = AA.box_eval(arb(s_lo), arb(b), arb(x), arb(y))
            v = r.get('rho', math.inf) if r.get('why') is None else math.inf
            if v > worst[0]: worst = (v, (b, x, y), r.get('why'))
    out.append(dict(s_lo=s_lo, L=L, Lneg=Lp, worst_rho=worst[0], at=worst[1], why=worst[2], falha=bool(worst[0] >= 1)))
    print('faixa [%g,%g) L*=%g:  com L*=%g  pior rho=%.3f em (beta, Re, Im)=%s %s -> %s' % (s_lo, s_hi, L, Lp, worst[0], worst[1], worst[2] or '',
          'FALHA (controle ok)' if worst[0] >= 1 else 'passa'), flush=True)
# (c) verificador completo (driver2) na primeira celula sigma de algumas faixas, com L* - 1: deve registrar falha
import driver2 as D
for s_lo, s_hi, L in [rows[0], rows[len(rows)//2], rows[-1]]:
    sb = min(s_hi, s_lo*1.1); C = D.Cert(s_lo, sb, L - 1, 100.0, maxev=60000); C.probe = False
    C.run(ratio=sb/s_lo*1.0000001)
    f = C.fails[0] if C.fails else None
    print('verificador, sigma in [%g,%g], L*=%g (tabela: %g): avaliacoes=%d, falhas=%d, primeira=%s' % (s_lo, sb, L - 1, L, C.nev, len(C.fails), f), flush=True)
    out.append(dict(driver=True, s_lo=s_lo, s_hi=sb, Lneg=L - 1, evals=C.nev, nfails=len(C.fails), first=str(f)))
for P0 in [100.0, 30.0]:
    w = math.nextafter(P0**-0.5, 1.0)
    r = AA.ann_eval(arb(288), arb(0), 0.0, w)
    print('cauda sigma=288, beta=0, P0=%g: ok=%s rho=%s why=%s' % (P0, r.get('ok'), r.get('rho'), r.get('why')))
    out.append(dict(tail_P0=P0, ok=r.get('ok'), rho=r.get('rho'), why=r.get('why')))
json.dump(out, open('negctl.json', 'w'), indent=1, default=str)
