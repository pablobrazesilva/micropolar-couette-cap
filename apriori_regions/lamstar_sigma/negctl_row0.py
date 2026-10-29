"""controle negativo do verificador completo na faixa [3.7498, 4.45) com L* = 13 (tabela: 14)"""
import json, driver2 as D
C = D.Cert(3.7498, 4.0, 13.0, 100.0, maxev=60000); C.probe = False
C.run(ratio=4.0/3.7498*1.0000001)
f = C.fails[0] if C.fails else None
print('verificador, sigma in [3.7498,4], L*=13 (tabela: 14): avaliacoes=%d, falhas=%d, primeira=%s' % (C.nev, len(C.fails), f), flush=True)
json.dump(dict(s_lo=3.7498, s_hi=4.0, Lneg=13, evals=C.nev, nfails=len(C.fails), first=str(f)), open('negctl_row0.json', 'w'))
