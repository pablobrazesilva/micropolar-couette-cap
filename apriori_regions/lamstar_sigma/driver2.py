"""Certificacao rigorosa (Arb) da regiao a priori, terna (1;1;0,5), faixa [s_lo, s_hi) com Lambda*:
   cobre {(sigma, beta, Lambda): sigma in [s_lo,s_hi], 0 <= beta <= alpha_max(sigma) sigma^{-1/3},
          Re Lambda >= -delta', Im Lambda <= S/2, e  (Re Lambda >= L*  ou  Im Lambda <= -L*)}
   (simetria Im Lambda -> S - Im Lambda cobre Im Lambda >= S/2) por
     (i) caixas (box_eval) com |Lambda_1| < P0 em algum ponto;  (ii) cauda |Lambda_1| >= P0 (ann_eval).
   Bissecao adaptativa: divide a direcao cuja contracao ao centro mais reduz rho (sondagem).
   uso: python3 driver2.py s_lo s_hi Lstar [P0] [tag] [ratio]"""
import sys, json, math, time
import ap_arb as AA, fa_table as FT
from flint import arb
DLT_F = math.nextafter(float(AA.DLT.upper()), 1.0)          # >= delta'
S1 = 8/3
def S_of(s): return 2*s**(1/3)
def rho_box(b):
    sa, sb, ba, bb, xa, xb, ya, yb = b
    return AA.box_eval(AA.iv(sa, sb), AA.iv(ba, bb), AA.iv(xa, xb), AA.iv(ya, yb))
def rho_tail(sa, sb, ba, bb, P_lo, P_hi):
    w_lo = 0.0 if P_hi == math.inf else P_hi**-0.5
    w_hi = math.nextafter(P_lo**-0.5, 1.0)
    return AA.ann_eval(AA.iv(sa, sb), AA.iv(ba, bb), w_lo, w_hi)
def split(b, k):
    b = list(b); a, c = b[2*k], b[2*k + 1]
    m = math.sqrt(a*c) if k == 0 else 0.5*(a + c)
    b1 = b[:]; b1[2*k + 1] = m; b2 = b[:]; b2[2*k] = m
    return tuple(b1), tuple(b2)
def collapse(b, k):
    b = list(b); a, c = b[2*k], b[2*k + 1]
    m = math.sqrt(a*c) if k == 0 else 0.5*(a + c); b[2*k] = b[2*k + 1] = m
    return tuple(b)
class Cert:
    def __init__(s, s_lo, s_hi, L, P0, maxev=3000000):
        s.s_lo, s.s_hi, s.L, s.P0 = s_lo, s_hi, L, P0
        s.gx = 0.4; s.db = 0.25; s.bg = 1.15; s.wsig = 0.05; s.wbet = 0.15; s.wre = 0.03; s.wim = 0.08; s.verbose = False; s.probe = True; s.nev = 0; s.ok = 0; s.fails = []; s.maxrho = 0.0; s.maxev = maxev; s.tail = 0; s.nbox = 0; s.skip = 0
    def far(s, b):
        """|Lambda_1| >= P0 em toda a caixa (coberta pela cauda)"""
        sa, sb, ba, bb, xa, xb, ya, yb = b
        re = xa + ba*ba + S1*sb**(-2/3)*(1 - 1e-12)
        im = 0.0 if ya <= 0 <= yb else min(abs(ya), abs(yb))
        return max(re, 0.0)**2 + im**2 >= s.P0**2*(1 + 1e-9)
    def widths(s, b):
        sa, sb, ba, bb, xa, xb, ya, yb = b
        dist = max(s.L, math.hypot(max(xa, 0.0), 0.0 if ya <= 0 <= yb else min(abs(ya), abs(yb))))
        return [math.log(sb/sa)/s.wsig, (bb - ba)*S_of(sa)/(s.wbet*max(1.0, 0.25*ba*S_of(sa))), (xb - xa)/(s.wre*dist), (yb - ya)/(s.wim*dist)]
    def run_box(s, b0):
        stack = [(b0, 0)]
        while stack:
            b, dep = stack.pop()
            if s.far(b): s.skip += 1; continue
            s.nev += 1
            if s.nev > s.maxev: s.fails.append(('limit', b)); return
            if s.verbose and s.nev % 500 == 0: print('    ev=%d ok=%d dep=%d box=%s' % (s.nev, s.nbox, dep, ['%.4g' % v for v in b]), flush=True)
            r = rho_box(b)
            if r.get('ok'):
                s.ok += 1; s.nbox += 1; s.maxrho = max(s.maxrho, r['rho']); continue
            if dep > 36: s.fails.append((r.get('why', 'rho'), r.get('rho'), b)); continue
            # sondagem: contrai cada direcao (nao-degenerada) ao centro
            best, kb = None, 0
            if not s.probe:
                w = s.widths(b); kb = max(range(4), key=lambda i: w[i]); best = 0
            for k in (range(4) if s.probe else []):
                if b[2*k] == b[2*k + 1]: continue
                s.nev += 1
                rk = rho_box(collapse(b, k)); v = rk.get('rho', 99.0) if 'rho' in rk else 99.0
                if best is None or v < best: best, kb = v, k
            c1, c2 = split(b, kb)
            stack.append((c1, dep + 1)); stack.append((c2, dep + 1))
    def run_tail(s, sa, sb, ba, bb):
        stack = [(sa, sb, ba, bb, s.P0, math.inf, 0)]
        while stack:
            sa, sb, ba, bb, Pl, Ph, dep = stack.pop()
            s.nev += 1; s.tail += 1
            r = rho_tail(sa, sb, ba, bb, Pl, Ph)
            if r.get('ok'): s.ok += 1; s.maxrho = max(s.maxrho, r['rho']); continue
            if dep > 30: s.fails.append(('tail ' + r.get('why', 'rho'), (sa, sb, ba, bb, Pl, Ph))); continue
            if (bb - ba)*S_of(sa) > 0.05:
                bm = 0.5*(ba + bb); stack += [(sa, sb, ba, bm, Pl, Ph, dep + 1), (sa, sb, bm, bb, Pl, Ph, dep + 1)]
            elif Ph == math.inf:
                stack += [(sa, sb, ba, bb, Pl, 4*Pl, dep + 1), (sa, sb, ba, bb, 4*Pl, math.inf, dep + 1)]
            else:
                sm = math.sqrt(sa*sb); stack += [(sa, sm, ba, bb, Pl, Ph, dep + 1), (sm, sb, ba, bb, Pl, Ph, dep + 1)]
    def grids(s, sa, sb):
        L, P0 = s.L, s.P0
        S = S_of(sa); bmax = math.nextafter(FT.alpha_max(sb)*sa**(-1/3), 10.0)
        bn = [0.0]; x = s.db/S
        while x < bmax:
            bn.append(x); x = x + s.db/S if x < 4/S else x*s.bg
        bn.append(bmax)
        Sh = S_of(sb)*(1 + 1e-12)
        xs = [L]; 
        while xs[-1] < P0: xs.append(xs[-1] + max(s.gx*xs[-1], 0.5))
        xs[-1] = P0
        ys = [Sh/2]; y = Sh/2
        while y > -P0:
            y = y - max(s.gx*max(abs(y), L), 0.5); ys.append(max(y, -P0))
        ys = ys[::-1]
        xI = [-DLT_F] + [L*k/4 for k in range(0, 5)]
        yI = [-L]; y = -L
        while y > -P0:
            y = y - max(s.gx*abs(y), 0.5); yI.append(max(y, -P0))
        yI = yI[::-1]
        return bn, xs, ys, xI, yI
    def run(s, ratio=1.03, log=None):
        n = max(1, math.ceil(math.log(s.s_hi/s.s_lo)/math.log(ratio)))
        edges = [s.s_lo*(s.s_hi/s.s_lo)**(k/n) for k in range(n + 1)]
        edges[0] = s.s_lo; edges[-1] = s.s_hi
        t0 = time.time()
        for i in range(n):
            sa, sb = edges[i], edges[i + 1]
            bn, xs, ys, xI, yI = s.grids(sa, sb)
            for j in range(len(bn) - 1):
                s.run_tail(sa, sb, bn[j], bn[j + 1])
            for j in range(len(bn) - 1):
                for a in range(len(xs) - 1):
                    for c in range(len(ys) - 1):
                        s.run_box((sa, sb, bn[j], bn[j + 1], xs[a], xs[a + 1], ys[c], ys[c + 1]))
                for a in range(len(xI) - 1):
                    for c in range(len(yI) - 1):
                        s.run_box((sa, sb, bn[j], bn[j + 1], xI[a], xI[a + 1], yI[c], yI[c + 1]))
                if s.fails: break
            print('  sigma in [%.4f, %.4f]: avaliacoes=%d, caixas ok=%d, cauda=%d, puladas=%d, falhas=%d, max rho=%.4f, %.0fs'
                  % (sa, sb, s.nev, s.nbox, s.tail, s.skip, len(s.fails), s.maxrho, time.time() - t0), flush=True)
            if s.fails: return
if __name__ == '__main__':
    s_lo, s_hi, L = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    P0 = float(sys.argv[4]) if len(sys.argv) > 4 else 100.0
    tag = sys.argv[5] if len(sys.argv) > 5 else '%g_%g_%g' % (s_lo, s_hi, L)
    ratio = float(sys.argv[6]) if len(sys.argv) > 6 else 1.03
    t = time.time(); C = Cert(s_lo, s_hi, L, P0); C.run(ratio)
    res = dict(s_lo=s_lo, s_hi=s_hi, Lstar=L, P0=P0, ratio=ratio, evals=C.nev, boxes_ok=C.nbox, tail_evals=C.tail, skipped=C.skip,
               fails=[str(f) for f in C.fails[:20]], nfails=len(C.fails), max_rho=C.maxrho, sec=time.time() - t)
    print('RESULTADO', json.dumps(res), flush=True)
    json.dump(res, open('cert_%s.json' % tag, 'w'))
