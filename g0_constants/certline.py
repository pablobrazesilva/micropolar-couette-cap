"""Certify a(eta), at(eta) lower bounds and |F| upper bounds on the line (b=0) with uniform boxes,
   plus tails; dumps per-box data for the P(S') bound.   usage: certline.py eta TL TR h"""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, json, time; sys.path.insert(0, _os.path.join(_R, 'g0_constants'))
from cl import *
from multiprocessing import Pool

def job(args):
    t1, t2, eta = args
    out = []; st = [(t1, t2)]
    while st:
        a1, a2 = st.pop()
        d = box_line(a1, a2, eta)
        if d is None:
            m = 0.5*(a1 + a2); st += [(a1, m), (m, a2)]; continue
        out.append((a1, a2, d['Flo'], d['Fup'], d['ulo'], d['utlo'], float(d['u0'].mid()), float(d['ut0'].mid())))
    return out

if __name__ == '__main__':
    eta = sys.argv[1]; TL = float(sys.argv[2]); TR = float(sys.argv[3]); h = float(sys.argv[4])
    T0 = time.time()
    n = int(round((TL + TR)/h)); edges = [-TL + h*k for k in range(n + 1)]
    with Pool(2) as P:
        res = P.map(job, [(edges[i], edges[i+1], eta) for i in range(n)], chunksize=16)
    boxes = sorted([b for r in res for b in r])
    # coverage check: boxes tile [-TL, TR] exactly (edges are dyadic multiples of h)
    assert boxes[0][0] == -TL and boxes[-1][1] == TR
    for b, c in zip(boxes, boxes[1:]): assert b[1] == c[0]
    tr = tail_side(TR, eta, 'right'); tl = tail_side(TL, eta, 'left')
    a_lo = min(min(b[4] for b in boxes), tr[0], tl[0])
    at_lo = min(min(b[5] for b in boxes), tr[1], tl[1])
    ib = min(range(len(boxes)), key=lambda i: boxes[i][4]); ibt = min(range(len(boxes)), key=lambda i: boxes[i][5])
    out = dict(eta=eta, TL=TL, TR=TR, h=h, nboxes=len(boxes), a_lo=a_lo, at_lo=at_lo,
               a_box=boxes[ib][:2], at_box=boxes[ibt][:2],
               a_centre_min=min(b[6] for b in boxes), at_centre_min=min(b[7] for b in boxes),
               Fmin_lo=min(b[2] for b in boxes), Fmax_up=max(b[3] for b in boxes),
               tail_right=tr, tail_left=tl, seconds=time.time() - T0)
    json.dump(dict(out, boxes=boxes), open(_wout('g0_constants', 'line_%s.json' % eta), 'w'))
    print(json.dumps(out, indent=1))
