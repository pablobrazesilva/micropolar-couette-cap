"""Rigorous verification of the coupling criterion on {sigma0 <= sigma <= sigma_h, 0 <= beta <= beta_A(sigma)}
by adaptive boxes in (log sigma, beta).  Pluggable: CRIT (Arb box criterion) and CHOOSE (float heuristic that picks
the free parameters lambda, th_rho, variant at a point) -- for crit3 replace both.
usage: driver.py N J G cA sigma0 sigma_h tag [coupling(1|0)]"""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, os, time, json, math
sys.path.insert(0, _os.path.join(_R, 'g1_wall_lemma')); sys.path.insert(0, _os.path.join(_R, 'g1_wall_lemma/heuristic'))
from multiprocessing import Pool
from arbcrit import crit_box, params
from arbprim import arb, iv
import crit2 as C2

def CHOOSE(NJG, s, b, coupling):
    """ordered candidate parameter choices: float-best overall, then float-best with the other lambda"""
    out = []
    sc, o = C2.best(*NJG, s, max(b, 1e-9), coupling=coupling)
    out.append(dict(lam=o['lam'], th_rho=o['th_rho'], cone=o['cone']))
    sc2, o2 = C2.best(*NJG, s, max(b, 1e-9), lams=(1.0 - o['lam'],), coupling=coupling)
    out.append(dict(lam=o2['lam'], th_rho=o2['th_rho'], cone=o2['cone']))
    return out

def CRIT(P, s1, s2, b1, b2, ch, coupling):
    return crit_box(P, iv(s1, s2), iv(b1, b2), ch['lam'], ch['th_rho'], ch['cone'], route2=(b2 <= 1.0), coupling=coupling)

# ---- plug-in interface (crit3): a module M may define
#   M.CRIT(P, s1, s2, b1, b2, choice, coupling) -> dict(margin=<rigorous lower bound over the box>, why=...)
#   M.cands(NJG, s, b, coupling)                -> iterable of 'choice' dicts (float heuristic, any order)
#   M.params(N, J, G)                           -> P   (optional; default arbcrit.params)
# selected with env var G1_CRIT=<module name> (default: this file's crit2 port).
_PLUG = None
if os.environ.get('G1_CRIT'):
    import importlib
    _PLUG = importlib.import_module(os.environ['G1_CRIT'])
    CRIT = _PLUG.CRIT
    if hasattr(_PLUG, 'params'): params = _PLUG.params

MAXEVAL = int(os.environ.get('G1_MAXEVAL', '600'))
BEDGES = [0, 0.002, 0.005, 0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 0.8, 1, 1.3, 1.7, 2.2, 3, 4, 5,
          6.5, 8.5, 11, 14, 18, 24, 32]

def column(args):
    NJG, cA, s1, s2, coupling = args
    P = params(*[repr(x) for x in NJG]) if False else params(*NJG_str(NJG))
    bA = s2**(1/6)/math.sqrt(cA)
    bmin = float(os.environ.get('G1_BMIN', '0'))          # v5: strip beta >= bmin only (K' reduction, G2 cluster v3)
    ed = [bmin] + [b for b in BEDGES if bmin < b < bA] + [bA*(1 + 1e-12)]
    stack = [(s1, s2, ed[i], ed[i+1]) for i in range(len(ed) - 1)]
    nbox = 0; neval = 0; fails = []; mins = (9, None)
    while stack:
        a1, a2, b1, b2 = stack.pop()
        if b1 > a2**(1/6)/math.sqrt(cA): continue          # box entirely above beta_A
        ok = False; tried = []
        def cands():
            if _PLUG is not None and hasattr(_PLUG, 'cands'):
                for (sp, bp) in ((math.sqrt(a1*a2), 0.5*(b1 + b2)), (a1, b1), (a1, b2)):
                    yield from _PLUG.cands(NJG, sp, bp, coupling)
                return
            for (sp, bp) in ((math.sqrt(a1*a2), 0.5*(b1 + b2)), (a1, b1), (a1, b2)):
                sc, o = C2.best(*NJG, sp, max(bp, 1e-9), coupling=coupling)
                yield dict(lam=o['lam'], th_rho=o['th_rho'], cone=o['cone'])
                sc2, o2 = C2.best(*NJG, sp, max(bp, 1e-9), lams=(1.0 - o['lam'],), coupling=coupling)
                yield dict(lam=o2['lam'], th_rho=o2['th_rho'], cone=o2['cone'])
        for ch in cands():
            if ch in tried: continue
            tried.append(ch)
            o = CRIT(P, a1, a2, b1, b2, ch, coupling); neval += 1
            if o['margin'] > 0:
                ok = True; nbox += 1
                if o['margin'] < mins[0]: mins = (o['margin'], (a1, a2, b1, b2, ch))
                break
        if ok: continue
        ws = math.log(a2/a1)/0.004; wb = (b2 - b1)/(0.002 + 0.1*b1)
        if math.log(a2/a1) < 1e-6 and b2 - b1 < 1e-7:
            fails.append((a1, a2, b1, b2, o.get('why'))); continue
        if len(fails) > 3 or neval > MAXEVAL:
            if neval > MAXEVAL: fails.append(('MAXEVAL', a1, a2, b1, b2))
            break
        if ws >= wb and math.log(a2/a1) >= 1e-6:
            m = math.sqrt(a1*a2); stack += [(a1, m, b1, b2), (m, a2, b1, b2)]
        else:
            m = 0.5*(b1 + b2); stack += [(a1, a2, b1, m), (a1, a2, m, b2)]
    return dict(col=(s1, s2), nbox=nbox, neval=neval, fails=fails, min=mins)

def NJG_str(NJG): return tuple(('%r' % x) for x in NJG)

if __name__ == '__main__':
    N, J, G = map(float, sys.argv[1:4]); cA = float(sys.argv[4]); s0 = float(sys.argv[5]); sh = float(sys.argv[6])
    tag = sys.argv[7]; coupling = int(sys.argv[8]) if len(sys.argv) > 8 else 1
    T0 = time.time()
    edges = [s0]
    while edges[-1] < sh:
        r = 1.02 if edges[-1] < 1.2*s0 else 1.05
        edges.append(min(sh, edges[-1]*r))
    jobs = [((N, J, G), cA, edges[i], edges[i+1], coupling) for i in range(len(edges) - 1)]
    with Pool(2) as Pp:
        res = Pp.map(column, jobs, chunksize=1)
    nb = sum(r['nbox'] for r in res); ne = sum(r['neval'] for r in res)
    fails = [f for r in res for f in r['fails']]
    mn = min(res, key=lambda r: r['min'][0])['min']
    out = dict(NJG=(N, J, G), cA=cA, sigma0=s0, sigma_h=sh, coupling=coupling, beta_min=float(os.environ.get('G1_BMIN', '0')), columns=len(jobs), boxes=nb, evals=ne,
               fails=fails, min_margin=mn[0], min_box=mn[1], seconds=time.time() - T0)
    json.dump(out, open(_wout('g1_wall_lemma', 'out_%s.json' % tag), 'w'), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))
