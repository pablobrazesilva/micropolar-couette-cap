"""controles negativos da continuacao (limites finais): para cada faixa com limite b_f > 0, desce beta (passo 0,05) a partir de
   b_f ate o criterio pontual na reta FALHAR (pior rho >= 1) em sigma = s_lo; registra o beta de falha.
   Mais: verificador completo (driver5) num retangulo abaixo do beta de falha (deve registrar falhas)."""
import sys, json
from scan_line4 import worst
import driver5 as D5
lims = json.load(open(sys.argv[1])); out = []; DRV = json.loads(sys.argv[3]) if len(sys.argv) > 3 else []
for s_lo, s_hi, b_f in lims:
    b = round(b_f - 0.05, 2); res = None
    while b >= -1e-9:
        w = worst(s_lo, max(b, 0.0))
        if w[0] >= 1: res = (max(b, 0.0), w[0], w[1]); break
        b = round(b - 0.05, 2)
    out.append(dict(band=[s_lo, s_hi], limite=b_f, beta_falha=res[0] if res else None, rho=res[1] if res else None, Im=res[2] if res else None))
    print('[%g,%g]: limite certificado %.2f; o criterio pontual na reta falha em beta=%s (rho=%s, Im=%s)' % (s_lo, s_hi, b_f, res and res[0], res and round(res[1], 3), res and res[2]), flush=True)
for s_lo, s_hi, bl, bh in (DRV if len(sys.argv) > 3 else []):
    C = D5.LineCert(s_lo, s_hi, bl, bh, maxev=30000); C.run(ratio=s_hi/s_lo*1.0000001, verbose=False)
    print('verificador [%g,%g]x[%g,%g]: avaliacoes=%d falhas=%d primeira=%s' % (s_lo, s_hi, bl, bh, C.nev, len(C.fails), str(C.fails[0])[:220] if C.fails else None), flush=True)
    out.append(dict(driver=[s_lo, s_hi, bl, bh], evals=C.nev, nfails=len(C.fails), first=str(C.fails[0]) if C.fails else None))
json.dump(out, open('negctl_cont2.json', 'w'), indent=1, default=str)
