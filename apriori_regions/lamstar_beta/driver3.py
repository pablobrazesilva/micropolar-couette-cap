"""Certificacao Arb da regiao a priori DEPENDENTE DE beta, terna (1;1;0,5): retangulo [s_lo,s_hi] x [b_lo,b_hi] com
   Lambda* = L (numero) ou L = None ('NONE': nenhum autovalor com Re Lambda >= -delta').
   Cobre {sigma, beta no retangulo, Re Lambda >= -delta', Im Lambda <= S/2, Lambda fora de B(L*)} (NONE: sem B) por
   caixas (box_eval: criterio (m_L, r_4) e, se beta > 0, (m_R, m_L); vale o melhor) em |Lambda_1| < P0 e cauda (ann_eval).
   uso: python3 driver3.py s_lo s_hi b_lo b_hi L|NONE out.json [ratio]"""
import sys, json, math, time, hashlib
import ap_arb as AA, driver2 as D2
DLT_F = D2.DLT_F; S_of = D2.S_of
class Cert(D2.Cert):
    def __init__(s, s_lo, s_hi, b_lo, b_hi, L, P0=None, maxev=3000000):
        if P0 is None: P0 = 100.0 if b_hi <= 3.2 else 400.0          # cauda: ann_eval exige |Lambda_1| >> (nu_1 - g) beta^2
        super().__init__(s_lo, s_hi, L if L is not None else 0.0, P0, maxev)
        s.b_lo, s.b_hi, s.none = b_lo, b_hi, L is None
    def widths(s, b):
        sa, sb, ba, bb, xa, xb, ya, yb = b
        dist = max(s.L, 1.0, math.hypot(max(xa, 0.0), 0.0 if ya <= 0 <= yb else min(abs(ya), abs(yb))))
        return [math.log(sb/sa)/s.wsig, (bb - ba)*S_of(sa)/(s.wbet*max(1.0, 0.25*ba*S_of(sa))), (xb - xa)/(s.wre*dist), (yb - ya)/(s.wim*dist)]
    def run_tail(s, sa, sb, ba, bb):
        stack = [(sa, sb, ba, bb, s.P0, math.inf, 0)]
        while stack:
            sa, sb, ba, bb, Pl, Ph, dep = stack.pop()
            s.nev += 1; s.tail += 1
            if s.nev > s.maxev: s.fails.append(('limit tail', (sa, sb, ba, bb))); return
            r = D2.rho_tail(sa, sb, ba, bb, Pl, Ph)
            if r.get('ok'): s.ok += 1; continue
            if dep > 24: s.fails.append(('tail ' + r.get('why', 'rho'), (sa, sb, ba, bb, Pl, Ph))); return
            if (bb - ba)*S_of(sa) > 0.05:
                bm = 0.5*(ba + bb); stack += [(sa, sb, ba, bm, Pl, Ph, dep + 1), (sa, sb, bm, bb, Pl, Ph, dep + 1)]
            elif Ph == math.inf:
                stack += [(sa, sb, ba, bb, Pl, 4*Pl, dep + 1), (sa, sb, ba, bb, 4*Pl, math.inf, dep + 1)]
            else:
                sm = math.sqrt(sa*sb); stack += [(sa, sm, ba, bb, Pl, Ph, dep + 1), (sm, sb, ba, bb, Pl, Ph, dep + 1)]
    def grids(s, sa, sb):
        L, P0 = s.L, s.P0
        S = S_of(sa)
        bn = [s.b_lo]; x = s.b_lo
        while True:
            x = x + s.db/S if x < 4/S else x*s.bg
            if x >= s.b_hi: break
            bn.append(x)
        bn.append(s.b_hi)
        Sh = S_of(sb)*(1 + 1e-12)
        if s.none:
            xs = [-DLT_F, 0.0, 0.25, 0.5, 1.0]
            while xs[-1] < P0: xs.append(xs[-1] + max(s.gx*xs[-1], 0.5))
            xs[-1] = P0
            ys = [Sh/2]; y = Sh/2
            while y > -P0:
                y = y - max(s.gx*abs(y), 0.5); ys.append(max(y, -P0))
            return bn, xs, ys[::-1], [], []
        xs = [L]
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
        return bn, xs, ys, xI, yI[::-1]
def md5(fn): return hashlib.md5(open(fn, 'rb').read()).hexdigest()
if __name__ == '__main__':
    s_lo, s_hi, b_lo, b_hi = map(float, sys.argv[1:5]); L = None if sys.argv[5] == 'NONE' else float(sys.argv[5])
    out = sys.argv[6]; ratio = float(sys.argv[7]) if len(sys.argv) > 7 else 1.1
    t = time.time(); C = Cert(s_lo, s_hi, b_lo, b_hi, L); C.probe = False; C.run(ratio)
    res = dict(NJG=[1, 1, 0.5], s_lo=s_lo, s_hi=s_hi, beta_lo=b_lo, beta_hi=b_hi, Lstar=('NONE' if L is None else L),
               delta_prime='0.05*0.5^(1/3)', P0=C.P0, ratio=ratio, evals=C.nev, boxes_ok=C.nbox, tail_evals=C.tail,
               skipped=C.skip, fails=[str(f) for f in C.fails[:20]], nfails=len(C.fails), sec=round(time.time() - t, 1),
               code={'ap_arb.py': md5('ap_arb.py'), 'driver2.py': md5('driver2.py'), 'driver3.py': md5('driver3.py')},
               ok=(len(C.fails) == 0))
    print('RESULTADO', json.dumps(res), flush=True)
    json.dump(res, open(out, 'w'), indent=1)
