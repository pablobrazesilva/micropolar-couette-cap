"""Driver: adaptive certification on the compact (t,rho) rectangle + tail A strip.
usage: python3 run_W.py eta th1 th2 mmin TL TR R tag
"""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, time, json, math
sys.path.insert(0, _os.path.join(_R, 'g0_sector_lemma'))
from multiprocessing import Pool
from g0lib import *

def certify_box(args):
    t1, t2, r1, r2, eta, th1, th2, mmin = args
    stack = [(t1, t2, r1, r2)]
    n = 0; fails = []; minm = [9, 9, 9]; neval = 0; minw = 9; ext = [9, -9, 9]
    while stack:
        a1, a2, b1, b2 = stack.pop()
        d = box_affine(a1, a2, b1, b2, eta); neval += 1
        ok, mg = check_box(d, th1, th2, mmin)
        if ok:
            n += 1
            g = lb(abs(d['G0']))
            minm[0] = min(minm[0], mg[0]/ub(abs(d['G0'])))
            minm[1] = min(minm[1], -mg[1]/ub(abs(d['G0'])))
            if mg[2] is not None: minm[2] = min(minm[2], mg[2])
            minw = min(minw, a2 - a1)
            ag = float((d['G0'].arg()/PI).mid()); ext = [min(ext[0], ag), max(ext[1], ag), min(ext[2], float(abs(d['G0']).mid()))]
            continue
        # point check at centre (detect genuine violation)
        if d is not None and d.get('k1') is not None:
            G0 = d['G0']
            v1 = (G0*(-I*PI*arb(th1)).exp()).imag; v2 = (G0*(-I*PI*arb(th2)).exp()).imag
            if ub(v1) < 0 or lb(v2) > 0 or (mmin is not None and ub(abs(G0)) < mmin):
                fails.append(('VIOLATION', a1, a2, b1, b2, str(G0)))
                continue
        if (a2 - a1) < 1e-6 and (b2 - b1) < 1e-6:
            fails.append(('TOO_SMALL', a1, a2, b1, b2)); continue
        st = d['st'] if d else 1; sr = d['sr'] if d else 1
        if (st >= sr and a2 - a1 > 1e-9) or b2 - b1 < 1e-9:
            m = 0.5*(a1 + a2); stack += [(a1, m, b1, b2), (m, a2, b1, b2)]
        else:
            m = 0.5*(b1 + b2); stack += [(a1, a2, b1, m), (a1, a2, m, b2)]
    return dict(box=(t1, t2, r1, r2), n=n, neval=neval, fails=fails, minm=minm, minw=minw, ext=ext)

def tailA_job(args):
    t1, t2, R, eta, th1, th2, mmin = args
    stack = [(t1, t2)]; n = 0; fails = []; worst = [9, -9, 1e300]
    while stack:
        a1, a2 = stack.pop()
        r = tailA_box(a1, a2, R, eta)
        if r is not None and r[0] > ub(arb(th1)) and r[1] < lb(arb(th2)) and (mmin is None or r[2] > mmin):
            n += 1; worst = [min(worst[0], r[0]), max(worst[1], r[1]), min(worst[2], r[2])]
        elif a2 - a1 > 1e-4:
            m = 0.5*(a1 + a2); stack += [(a1, m), (m, a2)]
        else:
            fails.append((a1, a2, r))
    return dict(n=n, fails=fails, worst=worst)

if __name__ == '__main__':
    eta = sys.argv[1]; th1 = sys.argv[2]; th2 = sys.argv[3]
    mmin = None if sys.argv[4] == 'none' else ub(arb(sys.argv[4]))  # rigorous: certified |G| > this >= exact decimal
    TL = float(sys.argv[5]); TR = float(sys.argv[6]); R = float(sys.argv[7]); tag = sys.argv[8]
    T0 = time.time()
    tB = tailB_check(eta, TL, TR, th1, th2, mmin)
    print('tailB:', tB, flush=True)
    tedges = [-TL + 0.5*k for k in range(int(round((TL + TR)/0.5)) + 1)]
    redges = [0, 0.05, 0.1, 0.2, 0.4, 0.7, 1, 1.5, 2, 3, 4, 6, 8, 11, 15, 20, 27, R]
    jobs = [(tedges[i], tedges[i+1], redges[j], redges[j+1], eta, th1, th2, mmin)
            for i in range(len(tedges)-1) for j in range(len(redges)-1)]
    jobsA = [(tedges[i], tedges[i+1], R, eta, th1, th2, mmin) for i in range(len(tedges)-1)]
    with Pool(2) as P:
        res = P.map(certify_box, jobs, chunksize=4)
        resA = P.map(tailA_job, jobsA)
    ntot = sum(r['n'] for r in res); nev = sum(r['neval'] for r in res)
    fails = [f for r in res for f in r['fails']]
    minm = [min(r['minm'][k] for r in res) for k in range(3)]
    nA = sum(r['n'] for r in resA); failsA = [f for r in resA for f in r['fails']]
    worstA = [min(r['worst'][0] for r in resA), max(r['worst'][1] for r in resA), min(r['worst'][2] for r in resA)]
    out = dict(eta=eta, th1=th1, th2=th2, mmin=mmin, TL=TL, TR=TR, R=R,
               compact_boxes=ntot, compact_evals=nev, compact_fails=fails,
               min_rel_margin_lower=minm[0], min_rel_margin_upper=minm[1], min_absG_lower=minm[2],
               centre_arg_min=min(r['ext'][0] for r in res), centre_arg_max=max(r['ext'][1] for r in res), centre_absG_min=min(r['ext'][2] for r in res),
               min_box_width_t=min(r['minw'] for r in res),
               tailA_boxes=nA, tailA_fails=failsA, tailA_arg_range=worstA[:2], tailA_absG_low=worstA[2],
               tailB=tB, seconds=time.time() - T0)
    json.dump(out, open(_wout('g0_sector_lemma', 'out_%s.json' % tag), 'w'), indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if 'fails' not in k}, indent=1, default=str))
    print('compact fails:', len(fails), fails[:10]); print('tailA fails:', len(failsA), failsA[:5])
