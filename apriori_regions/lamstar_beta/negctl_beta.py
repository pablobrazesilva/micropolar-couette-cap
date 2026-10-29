"""controles negativos da tabela Lambda*(sigma, beta).
 (a) retangulos com Lambda*: com Lambda* - 0.5, pior rho pontual (sigma = s_lo, beta em {b_lo, meio, b_hi}) na nova fronteira;
 (b) retangulos NONE: com beta = b_lo - 0.1 (sigma = s_lo e s_hi), pior rho pontual na regiao inteira (sem B);
 (c) verificador completo (driver3) num retangulo NONE com b_lo - 0.15: deve falhar.
 uso: python3 negctl_beta.py table.json out.json"""
import sys, json, math, ap_arb as AA
from flint import arb
from scan_pts import pts
from scan_vazio import worst as worst_full
import driver3 as D3
rects = json.load(open(sys.argv[1])); out = []
def worst_L(sg, bs, Ls):
    S = 2*sg**(1/3); w = (0, None)
    for b in bs:
        for L in pts(Ls, S):
            o = AA.box_eval(arb(sg), arb(b), arb(max(L.real, -0.0396)), arb(L.imag))
            v = o.get('rho', math.inf) if o.get('why') is None else math.inf
            if v > w[0]: w = (v, (b, max(L.real, -0.0396), L.imag))
    return w
for s_lo, s_hi, b_lo, b_hi, L in rects:
    if L == 'NONE':
        b = b_lo - 0.1
        w1 = worst_full(s_lo, b); w2 = worst_full(s_hi, b)
        w = max(w1[0], w2[0])
        out.append(dict(rect=[s_lo, s_hi, b_lo, b_hi, L], test='beta=b_lo-0.1', worst_rho=w, at=str(w1[1] if w1[0] >= w2[0] else w2[1]), falha=w >= 1))
        print('NONE [%g,%g]x[%g,%.3f]: com beta=%.2f pior rho=%s -> %s' % (s_lo, s_hi, b_lo, b_hi, b, w if w < 9 else 'falha(E/q)', 'FALHA (ok)' if w >= 1 else 'passa'), flush=True)
    else:
        res = {}
        for dL in (0.5, 1.0):
            if L - dL <= 0:                                       # sem caixa B: regiao inteira (como NONE) em beta = b_lo
                w = worst_full(s_lo, b_lo); res[dL] = (w[0] if w[0] < 9 else math.inf, ('B vazio', str(w[1]))); continue
            res[dL] = worst_L(s_lo, [b_lo, 0.5*(b_lo + b_hi), b_hi], L - dL)
        out.append(dict(rect=[s_lo, s_hi, b_lo, b_hi, L], test='Lstar-0.5/-1', worst_rho_m05=res[0.5][0], worst_rho_m1=res[1.0][0],
                        at_m1=str(res[1.0][1]), falha_m1=res[1.0][0] >= 1))
        print('[%g,%g]x[%g,%g] L*=%g: com L*-0.5 pior rho=%.3f; com L*-1 pior rho=%.3f -> %s' % (s_lo, s_hi, b_lo, b_hi, L, res[0.5][0], res[1.0][0],
              'FALHA (ok)' if res[1.0][0] >= 1 else 'passa'), flush=True)
# (c)
for s_lo, s_hi, b_lo, b_hi, L in [r for r in rects if r[4] == 'NONE'][::7]:
    C = D3.Cert(s_lo, min(s_hi, s_lo*1.1), b_lo - 0.15, b_lo, None, maxev=40000); C.probe = False
    C.run(min(s_hi, s_lo*1.1)/s_lo*1.0000001)
    print('verificador NONE em [%g,%g]x[%g,%g]: avaliacoes=%d falhas=%d primeira=%s' % (s_lo, min(s_hi, s_lo*1.1), b_lo - 0.15, b_lo, C.nev, len(C.fails), str(C.fails[0])[:200] if C.fails else None), flush=True)
    out.append(dict(rect=[s_lo, s_hi, b_lo - 0.15, b_lo, 'NONE'], test='driver', evals=C.nev, nfails=len(C.fails), first=str(C.fails[0]) if C.fails else None))
json.dump(out, open(sys.argv[2], 'w'), indent=1, default=str)
