"""controle negativo do verificador (driver5) com orcamento pequeno: no retangulo abaixo do limite, o verificador
   NAO certifica (falha por caixa ou por orcamento)."""
import sys, json, driver5 as D5
s_lo, s_hi, bl, bh, mx = map(float, sys.argv[1:6])
C = D5.LineCert(s_lo, s_hi, bl, bh, maxev=int(mx)); C.run(ratio=s_hi/s_lo*1.0000001, verbose=False)
print('verificador [%g,%g]x[%g,%g] (orcamento %d): avaliacoes=%d, falhas=%d, primeira=%s' % (s_lo, s_hi, bl, bh, mx, C.nev, len(C.fails), str(C.fails[0])[:250] if C.fails else None), flush=True)
