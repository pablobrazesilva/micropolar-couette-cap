"""Certificacao 'NONE por continuacao' (Teorema de continuacao, relatorio apriori_beta_v2 §3), terna (1;1;0,5).

   Retangulo Q = [s_lo, s_hi] x [b_lo, b_hi] adjacente (em beta = b_hi) a um retangulo NONE ja certificado.
   Verifica-se, para TODO (sigma, beta) em Q:
     (a) nenhum autovalor na RETA Re Lambda = -delta', Im Lambda in [-P0, S/2]  (a simetria Im -> S - Im cobre o resto),
         por caixas (sigma, beta, Im Lambda) com Re Lambda na bola de -delta': criterio LG (ap_arb.box_eval) ou,
         se este falhar, o criterio Airy exato (ap_airy3.airy3_box);
     (b) a cauda |Lambda_1| >= P0, Re Lambda >= -delta' (ap_arb.ann_eval).
   Entao o numero de autovalores em {Re Lambda > -delta'} e' constante em Q (continuidade da funcao de Evans, Hurwitz)
   e igual ao do retangulo NONE adjacente: zero.
   uso: python3 driver5.py s_lo s_hi b_lo b_hi out.json"""
import sys, json, math, time, hashlib
import ap_arb as AA, ap_airy3 as A3, driver2 as D2
from flint import arb
DLTB = -AA.DLT                                       # bola que contem -delta' exato
def S_of(s): return 2*s**(1/3)
class LineCert:
    def __init__(s, s_lo, s_hi, b_lo, b_hi, P0=None, maxev=2000000):
        s.s_lo, s.s_hi, s.b_lo, s.b_hi = s_lo, s_hi, b_lo, b_hi
        s.P0 = P0 if P0 is not None else (100.0 if b_hi <= 3.2 else 400.0)
        s.nev = 0; s.nLG = 0; s.nAE = 0; s.fails = []; s.tail = 0; s.maxev = maxev; s.ok = 0
    def ev(s, b):
        sa, sb, ba, bb, ya, yb = b
        sig, bet, Y = AA.iv(sa, sb), AA.iv(ba, bb), AA.iv(ya, yb)
        r = AA.box_eval(sig, bet, DLTB, Y)
        if r.get('ok'): s.nLG += 1; return r
        best = r
        for rows in (('mR', 'mL'), ('mL', 'r4')):           # Airy exato: linhas (m_R, m_L) [beta > 0] ou (m_L, r_4')
            if rows == ('mR', 'mL') and not ba > 0: continue
            r2 = A3.airy3_safe(sig, bet, DLTB, Y, rows=rows); s.nAE += 1
            if r2.get('ok'): return r2
            if 'rho' in r2 and r2['rho'] < best.get('rho', math.inf): best = r2
        return best
    def widths(s, b):
        sa, sb, ba, bb, ya, yb = b
        return [math.log(sb/sa)/0.05, (bb - ba)/0.03, (yb - ya)/(0.15*max(1.0, min(abs(ya), abs(yb))/3))]
    def run_box(s, b0):
        stack = [(b0, 0)]
        while stack:
            b, dep = stack.pop()
            s.nev += 1
            if s.nev > s.maxev: s.fails.append(('limit', b)); return
            r = s.ev(b)
            if r.get('ok'): s.ok += 1; continue
            if dep > 30: s.fails.append((r.get('why', 'rho'), r.get('rho'), b)); continue
            w = s.widths(b); k = max(range(3), key=lambda i: w[i])
            b = list(b); a, c_ = b[2*k], b[2*k + 1]; m = math.sqrt(a*c_) if k == 0 else 0.5*(a + c_)
            b1 = b[:]; b1[2*k + 1] = m; b2 = b[:]; b2[2*k] = m
            stack.append((tuple(b1), dep + 1)); stack.append((tuple(b2), dep + 1))
    def run_tail(s, sa, sb, ba, bb):
        stack = [(sa, sb, ba, bb, s.P0, math.inf, 0)]
        while stack:
            sa, sb, ba, bb, Pl, Ph, dep = stack.pop()
            s.nev += 1; s.tail += 1
            r = D2.rho_tail(sa, sb, ba, bb, Pl, Ph)
            if r.get('ok'): continue
            if dep > 24: s.fails.append(('tail ' + r.get('why', 'rho'), (sa, sb, ba, bb, Pl, Ph))); return
            if (bb - ba)*S_of(sa) > 0.05:
                bm = 0.5*(ba + bb); stack += [(sa, sb, ba, bm, Pl, Ph, dep + 1), (sa, sb, bm, bb, Pl, Ph, dep + 1)]
            elif Ph == math.inf:
                stack += [(sa, sb, ba, bb, Pl, 4*Pl, dep + 1), (sa, sb, ba, bb, 4*Pl, math.inf, dep + 1)]
            else:
                sm = math.sqrt(sa*sb); stack += [(sa, sm, ba, bb, Pl, Ph, dep + 1), (sm, sb, ba, bb, Pl, Ph, dep + 1)]
    def run(s, ratio=1.1, verbose=True):
        n = max(1, math.ceil(math.log(s.s_hi/s.s_lo)/math.log(ratio)))
        edges = [s.s_lo*(s.s_hi/s.s_lo)**(k/n) for k in range(n + 1)]; edges[0] = s.s_lo; edges[-1] = s.s_hi
        nb = max(1, math.ceil((s.b_hi - s.b_lo)/0.05))
        bn = [s.b_lo + (s.b_hi - s.b_lo)*k/nb for k in range(nb + 1)]; bn[0] = s.b_lo; bn[-1] = s.b_hi
        t0 = time.time()
        for i in range(n):
            sa, sb = edges[i], edges[i + 1]
            Sh = S_of(sb)*(1 + 1e-12)
            ys = [Sh/2]; y = Sh/2
            while y > -s.P0:
                y = y - max(0.25, 0.25*abs(y)); ys.append(max(y, -s.P0))
            ys = ys[::-1]
            for j in range(nb):
                s.run_tail(sa, sb, bn[j], bn[j + 1])
                for a in range(len(ys) - 1):
                    s.run_box((sa, sb, bn[j], bn[j + 1], ys[a], ys[a + 1]))
                if s.fails: break
            if verbose:
                print('  sigma in [%.4f, %.4f]: avaliacoes=%d (LG ok %d, Airy %d), cauda=%d, falhas=%d, %.0fs'
                      % (sa, sb, s.nev, s.nLG, s.nAE, s.tail, len(s.fails), time.time() - t0), flush=True)
            if s.fails: return
def md5(fn): return hashlib.md5(open(fn, 'rb').read()).hexdigest()
if __name__ == '__main__':
    s_lo, s_hi, b_lo, b_hi = map(float, sys.argv[1:5]); out = sys.argv[5]
    t = time.time(); C = LineCert(s_lo, s_hi, b_lo, b_hi); C.run()
    res = dict(NJG=[1, 1, 0.5], s_lo=s_lo, s_hi=s_hi, beta_lo=b_lo, beta_hi=b_hi, Lstar='NONE', metodo='continuacao (reta + cauda)',
               delta_prime='0.05*0.5^(1/3)', P0=C.P0, evals=C.nev, LG_ok=C.nLG, airy_evals=C.nAE, tail_evals=C.tail,
               fails=[str(f) for f in C.fails[:10]], nfails=len(C.fails), sec=round(time.time() - t, 1),
               code={f: md5(f) for f in ('ap_arb.py', 'ap_airy.py', 'ap_airy3.py', 'driver2.py', 'driver5.py')}, ok=(len(C.fails) == 0))
    print('RESULTADO', json.dumps(res), flush=True)
    json.dump(res, open(out, 'w'), indent=1)
