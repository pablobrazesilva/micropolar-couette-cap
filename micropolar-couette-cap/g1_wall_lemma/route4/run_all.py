"""Route 4: certify the joint constants on all beta nodes and assemble the certificate (resumable).
  python3 run_all.py            # builds the Arb tables in $CAP_WORK/route4 if needed, runs u/v x exp/s, writes M_cert_1_1_0.5.json
  python3 run_all.py u exp      # one family only
The threshold passed to mjoint_cert.py only decides which rows get the 2-D treatment (any value is valid)."""
import sys, json, os, subprocess, hashlib, time
HERE = os.path.dirname(os.path.abspath(__file__)); R0 = os.path.abspath(os.path.join(HERE, '..', '..'))
W = os.path.join(os.environ.get('CAP_WORK', os.path.join(R0, 'work')), 'route4'); os.makedirs(W, exist_ok=True)
NODES = [0, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16,
         0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
LINE, FTAB = os.path.join(W, 'line_tab.npz'), os.path.join(W, 'Ftab128.npz')
FT = json.load(open(os.path.join(R0, 'g1_wall_lemma', 'heuristic', 'Mtab_1_1_0.5.json')))     # floating-point values (thresholds only)
fb = sorted((float(k), v) for k, v in FT.items() if k != 's')
def guess(pair, wk, b):
    i = 0 if pair == 'u' else 1
    if wk == 's': return FT['s'][i]
    lo = max(x for x in fb if x[0] <= b); hi = min(x for x in fb if x[0] >= b)
    return lo[1][i] if hi[0] == lo[0] else lo[1][i] + (hi[1][i] - lo[1][i])*(b - lo[0])/(hi[0] - lo[0])
def run(pair, wk):
    out = os.path.join(W, 'M_%s_%s.jsonl' % (pair, wk)); done = set()
    if os.path.exists(out): done = {json.loads(l)['beta'] for l in open(out)}
    for b in (NODES if wk == 'exp' else [0]):
        if b in done: continue
        p = subprocess.run([sys.executable, os.path.join(HERE, 'mjoint_cert.py'), pair, wk, repr(b), repr(0.96*guess(pair, wk, b)/1.02),
                            LINE, FTAB], capture_output=True, text=True)
        line = [l for l in p.stdout.splitlines() if l.startswith('RESULT ')]
        if not line: print('FAILED', pair, wk, b, p.stderr[-800:], flush=True); continue
        r = json.loads(line[0][7:])
        with open(out, 'a') as f: f.write(json.dumps(r) + '\n')
        print(pair, wk, b, 'M <= %.5f (lower %.5f) %.0fs' % (r['M'], r['M_lower'], r['seconds']), flush=True)
def md5(p): return hashlib.md5(open(p, 'rb').read()).hexdigest()
def sha(p):
    """sha256 of the array contents of an .npz (independent of the zip container metadata)"""
    import numpy as np
    d = np.load(p); h = hashlib.sha256()
    for k in sorted(d.files): h.update(k.encode()); h.update(np.ascontiguousarray(d[k]).tobytes())
    return h.hexdigest()
def roundup(x, n=6):
    e = 10**(n - 1 - int(__import__('math').floor(__import__('math').log10(abs(x)))))
    return __import__('math').ceil(x*e)/e
def assemble(dest):
    rec = {}
    for pair in ('u', 'v'):
        for wk in ('exp', 's'):
            for l in open(os.path.join(W, 'M_%s_%s.jsonl' % (pair, wk))):
                r = json.loads(l); rec[(pair, wk, r['beta'])] = r
    ok = all((p, 'exp', b) in rec for p in 'uv' for b in NODES) and all((p, 's', 0) in rec for p in 'uv')
    for r in rec.values():
        ok = ok and r['M'] >= max(r['regions'].values()) and r['M_lower'] <= r['M']
    M_exp = {repr(b): [roundup(rec[('u', 'exp', b)]['M']), roundup(rec[('v', 'exp', b)]['M'])] for b in NODES}
    lower = {repr(b): [rec[('u', 'exp', b)]['M_lower'], rec[('v', 'exp', b)]['M_lower']] for b in NODES}
    d = dict(NJG=[1, 1, 0.5], eta=0.05, definition='M = sup_{(t_k,t_f) in R^2} |f_beta|, Lemma F5; u: kernel g, forcing nu1; v: kernel nu1, forcing g',
             M_exp=M_exp, M_s=[roundup(rec[('u', 's', 0)]['M']), roundup(rec[('v', 's', 0)]['M'])],
             M_lower_exp=lower, M_lower_s=[rec[('u', 's', 0)]['M_lower'], rec[('v', 's', 0)]['M_lower']],
             runs={'%s_%s_%r' % k: dict(M=v['M'], M_lower=v['M_lower'], regions=v['regions'], strip_tk=v['strip_tk'], seconds=v['seconds'])
                   for k, v in sorted(rec.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2]))},
             code={f: md5(os.path.join(HERE, f)) for f in ('mjoint_cert.py', 'linetab.py', 'ftab.py', 'run_all.py')},
             tables={'line_tab.npz': sha(LINE), 'Ftab128.npz': sha(FTAB)}, all_ok=bool(ok))
    json.dump(d, open(dest, 'w'), indent=1)
    print('assembled', dest, 'all_ok', ok)
if __name__ == '__main__':
    if not os.path.exists(LINE): subprocess.run([sys.executable, os.path.join(HERE, 'linetab.py'), '12800', LINE], check=True)
    if not os.path.exists(FTAB): subprocess.run([sys.executable, os.path.join(HERE, 'ftab.py'), '-2560', '3840', FTAB], check=True)
    fams = [tuple(sys.argv[1:3])] if len(sys.argv) > 2 else [('u', 'exp'), ('v', 'exp'), ('u', 's'), ('v', 's')]
    for fam in fams: run(*fam)
    if len(sys.argv) <= 2: assemble(os.path.join(W, 'M_cert_1_1_0.5.json'))
