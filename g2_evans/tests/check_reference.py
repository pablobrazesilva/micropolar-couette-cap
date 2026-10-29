"""Valida o binario C contra valores de referencia da versao Python (python-flint, ArbEvansC, 256 bits, grau 40, h|A|~1):
   para cada ponto, |D_C - D_py| <= rad_C + rad_py + 1e-20 |D| (as duas sao bolas rigorosas: devem se intersectar)."""
import json, subprocess, os, re, sys
here = os.path.dirname(os.path.abspath(__file__)); exe = os.path.join(here, '..', 'evans')
ref = json.load(open(os.path.join(here, 'ref_python.json')))
bad = 0
from decimal import Decimal, getcontext
getcontext().prec = 50
for r in ref:
    out = subprocess.run([exe, 'point', '1', '1', '0.5', repr(float(r['alpha'])), repr(float(r['sigma'])), repr(float(r['Lam'][0])), repr(float(r['Lam'][1]))],
                         capture_output=True, text=True).stdout
    m = re.search(r'D_0 = (\S+) ([+-]) (\S+)\*I\s+rad\(re,im\) = (\S+) (\S+)', out)
    re_c = Decimal(m.group(1)); im_c = Decimal(m.group(3)) * (1 if m.group(2) == '+' else -1)
    radc = float(m.group(4)) + float(m.group(5))
    t = float(re.search(r'time=(\S+)', out).group(1))
    dre = abs(re_c - Decimal(r['re'])); dim = abs(im_c - Decimal(r['im']))
    mag = abs(Decimal(r['re'])) + abs(Decimal(r['im']))
    tol = Decimal(radc + r['rad']) + mag * Decimal('1e-18')    # C imprime 20 digitos: tolerancia de impressao
    ok = dre <= tol and dim <= tol
    bad += not ok
    print('alpha=%-5g sigma=%-5g Lam=%s: rel.dif=%.2e  rad_C=%.1e rad_py=%.1e  t_C=%.3fs t_py=%.3fs  %s'
          % (r['alpha'], r['sigma'], r['Lam'], float((dre + dim)/mag), radc, r['rad'], t, r['t_python'], 'OK' if ok else 'FALHA'))
sys.exit(1 if bad else 0)
