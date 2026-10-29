"""roda driver2 em sequencia para uma lista de faixas 's_lo:s_hi:Lstar' (celulas sigma de razao <= ratio).
   grava cert_<s_lo>_<s_hi>_<L>.json por faixa e um log por celula sigma."""
import sys, json, time, math
import driver2 as D
ratio = float(sys.argv[1]); out = sys.argv[2]
for spec in sys.argv[3:]:
    s_lo, s_hi, L = map(float, spec.split(':'))
    t = time.time(); C = D.Cert(s_lo, s_hi, L, 100.0); C.probe = False
    C.run(ratio)
    res = dict(s_lo=s_lo, s_hi=s_hi, Lstar=L, P0=100.0, ratio=ratio, evals=C.nev, boxes_ok=C.nbox, tail_evals=C.tail,
               skipped=C.skip, fails=[str(f) for f in C.fails[:20]], nfails=len(C.fails), max_rho=C.maxrho, sec=time.time() - t)
    print('RESULTADO', json.dumps(res), flush=True)
    with open(out, 'a') as f: f.write(json.dumps(res) + '\n')
