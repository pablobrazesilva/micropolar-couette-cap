#!/usr/bin/env python3
"""Level-1 check (Python standard library only, a few seconds): reread every recorded certificate on the logical path
of the Main Theorem and check the claims the article makes about it (no failures, triple (1,1,1/2), ranges that cover
what the proof needs, values used <= values certified, recorded code hashes = the code in this repository, and the
SHA-256 values quoted in the article).  It does NOT redo the interval computations (see scripts/recompute_samples.py).
The Evans-function sweep (G2) is checked separately by `g2_evans/sweep.py verify` (scripts/quick_verify.sh)."""
import glob, hashlib, json, math, os, sys

R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
C = os.path.join(R, 'certificates')
NFAIL = 0

def check(label, cond, detail=''):
    global NFAIL
    if not cond: NFAIL += 1
    print('  [%s] %s%s' % ('ok' if cond else 'FAIL', label, (' -- ' + detail) if detail else ''))

def J(*p): return json.load(open(os.path.join(C, *p)))
def sha256(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
def md5(p): return hashlib.md5(open(p, 'rb').read()).hexdigest()

print('Article hashes (Remark G2cert)')
for p, h in (('g2_evans/verify_final.txt', 'ba84b458863ccc7a3f2afdd5d62d4e435e6d0a5a492f36dc420c741c9d9f32a3'),
             ('apriori_regions/lamstar.txt', '77907aa1deadf37a'), ('apriori_regions/lamstar_beta.txt', 'b8a3e384f93bb87f')):
    check('sha256(certificates/%s) = %s...' % (p, h[:16]), sha256(os.path.join(C, p)).startswith(h))
man = [l for l in open(os.path.join(C, 'g2_evans/ARCHIVE_MANIFEST.txt'))]
arch = [l.split()[-1] for l in man if l.startswith('# sha256 of the archive')]
mem = [l.split() for l in man if l.strip() and not l.startswith('#')]
here = [m for m in mem if m[4] != 'NOT']
badm = [m[2] for m in here if sha256(os.path.join(R, m[4])) != m[0]]
check('archive g2_certificados_final.tar.gz (sha256 0b3fde2f...): %d of its %d files are in this repository, byte for byte'
      % (len(here), len(mem)), arch == ['0b3fde2f7cb27b4773e1ea7466ee8ea1a8f3096a71232b02548b56b80429e4aa'] and not badm, ', '.join(badm[:3]))
for p, h in (('g2_evans/sweep.py', '920145f9'), ('g2_evans/src/evans.c', '176d5693')):
    check('sha256(%s) = %s...' % (p, h), sha256(os.path.join(R, p)).startswith(h))

print('Copies used as inputs by g2_evans/ are identical to certificates/')
pairs = [('g2_evans/lamstar.txt', 'certificates/apriori_regions/lamstar.txt'),
         ('g2_evans/lamstar_beta.txt', 'certificates/apriori_regions/lamstar_beta.txt'),
         ('g2_evans/fa_cert.json', 'certificates/apriori_regions/fa_cert.json'),
         ('g2_evans/g1strip.txt', 'certificates/g1_wall_lemma/strips/g1strip.txt')]
pairs += [('g2_evans/apriori_cert/' + os.path.basename(f), 'certificates/apriori_regions/apriori_cert/' + os.path.basename(f))
          for f in glob.glob(os.path.join(C, 'apriori_regions/apriori_cert/*.json'))]
pairs += [('g2_evans/g1cert/' + os.path.basename(f), 'certificates/g1_wall_lemma/strips/' + os.path.basename(f))
          for f in glob.glob(os.path.join(C, 'g1_wall_lemma/strips/out_S_*.json'))]
bad = [a for a, b in pairs if not (os.path.exists(os.path.join(R, a)) and md5(os.path.join(R, a)) == md5(os.path.join(R, b)))]
check('%d input files identical' % len(pairs), not bad, ', '.join(bad[:5]))

print('Appendix A: sector lemma (cone) at eta = 0.05 and 0.1, all rho >= 0')
for f, eta, th1, th2, m in (('out_C1_005.json', '0.05', -0.0836, 0.54937, 0.6083), ('out_C1_010.json', '0.1', -0.0837, 0.56137, 0.5827)):
    d = J('g0_sector_lemma', f)
    ok = (d['eta'] == eta and float(d['th1']) == th1 and float(d['th2']) == th2 and abs(float(d['mmin']) - m) < 1e-12
          and d['compact_fails'] == [] and d['tailA_fails'] == [] and all(x[-1] is True for x in d['tailB'])
          and float(d['TL']) == 20 and float(d['R']) == 40)
    check('eta=%s: arg G in (%g pi, %g pi), |G| > %g; compact %d boxes, tail A %d boxes, tail B both sides'
          % (eta, th1, th2, m, d['compact_boxes'], d['tailA_boxes']), ok)
    if eta == '0.05':
        c = math.cos((th2 - th1)*math.pi/2)
        check('cos(theta) = cos((th2-th1) pi/2) = %.7f >= 0.545117' % c, c >= 0.545117)

print('Appendix B/D: boundary constants and tables (eta = 0.05)')
d = J('g0_constants', 'line_0.05.json')
xs = sorted((b[0], b[1]) for b in d['boxes'])
contig = xs[0][0] == -20 and xs[-1][1] == 30 and all(xs[i][1] == xs[i + 1][0] for i in range(len(xs) - 1))
check('line_0.05.json: %d contiguous boxes covering [-20,30], tails both sides' % len(xs), contig and d['TL'] == 20 and d['TR'] == 30)
P = [(float(a), float(b)) for a, b in J('g0_constants', 'ptable_0.05.json')]
check("P(S') table: %d nodes on [%g,%g], non-increasing" % (len(P), P[0][0], P[-1][0]),
      all(P[i][1] >= P[i + 1][1] for i in range(len(P) - 1)) and P[0][0] == 1 and P[-1][0] == 26)
t12, t36 = J('g0_constants', 'tp_cert_0.05.json'), J('g0_constants', 'T36_cert_0.05.json')
for nm, d in (('OmA', t12), ('OmPsi', t12), ('KPi', t36), ('KPi1', t36), ('Q', t36), ('Qs', t36)):
    v = d[nm]
    check('table %-5s: %d nodes on l in [0,16], non-increasing majorant' % (nm, len(v)),
          len(v) == 641 and all(v[i] >= v[i + 1] for i in range(len(v) - 1)))
for f in ('T8_cert_1_1_0.5.json', 'T8_cert_1_1_0.5_s215.json'):
    d = J('g0_constants', f)
    fin = [e for e in d['Theta_edge'] if e[2] < 10]
    check('%s: sigma_min=%g, %d of %d beta-cells certified' % (f, d['sigma_min'], len(fin), len(d['Theta_edge'])), len(fin) > 0)
d = J('g0_constants', 'J_cert_1_1_0.5.json')
check('J_cert_1_1_0.5.json: route-3 constants for (g,nu1) and (nu1,g)', set(d) >= {'g1', '1g'})

print('Appendix D: Pi constants with tails in Arb')
d = J('g1_wall_lemma', 'cert_pi_arb.json'); p0, p1 = J('g1_wall_lemma', 'cert_pi_0.05.json'), J('g1_wall_lemma', 'cert_pipsi.json')
used = dict(Pi=p0['Pi'], Pi1=p0['Pi1'], Pi_psi=p1['Pi_psi'], Pi1_psi=p1['Pi1_psi'])
check('values used by the G1 criterion = values listed in cert_pi_arb.json', used == d['used_in_G1'])
check('every value used >= its new Arb certificate (tails included)',
      d['all_used_are_upper_bounds'] is True and all(d['new'][k] <= used[k] for k in used),
      ', '.join('%s %.8f <= %.8f' % (k, d['new'][k], used[k]) for k in used))

print('Appendix D: Theorem C criterion K > 0 for (1,1,1/2), sigma >= 288')
runs = [('out_C4R2b.json', 288, 300), ('out_C4R2.json', 300, 650), ('out_C3R2.json', 650, 1100), ('out_C3A2.json', 1100, 1900),
        ('out_A2.json', 1900, 3e5)]
for f, lo, hi in runs:
    d = J('g1_wall_lemma', f)
    ok = (d['fails'] == [] and [float(x) for x in d['NJG']] == [1, 1, 0.5] and d['cA'] == 0.796 and d['coupling'] == 1
          and d['sigma0'] == lo and d['sigma_h'] == hi and d.get('beta_min', 0) == 0)
    check('[%g, %g]: %d boxes, 0 failures, min margin %.3g' % (lo, hi, d['boxes'], d['min_margin']), ok)
check('the five ranges are contiguous and cover [288, 3e5]; sigma >= 3e5 by the tail lemma (tail.py)',
      all(runs[i][2] == runs[i + 1][1] for i in range(len(runs) - 1)))
strips = [l.split('#')[0].split() for l in open(os.path.join(C, 'g1_wall_lemma/strips/g1strip.txt'))]
for lo, hi, bmin, fn in [s for s in strips if len(s) >= 4]:
    d = J('g1_wall_lemma', 'strips', os.path.basename(fn))
    ok = (d['fails'] == [] and d['sigma0'] <= float(lo) and d['sigma_h'] >= float(hi) and d['beta_min'] <= float(bmin)
          and d['cA'] == 0.796 and d['coupling'] == 1)
    check('band sigma in [%s,%s], beta >= %s: %d boxes, 0 failures' % (lo, hi, bmin, d['boxes']), ok)

print('Remark cA and Appendix E: F_A, Lambda*(sigma), Lambda*(sigma,beta)')
d = J('apriori_regions', 'fa_cert.json')
pairs = sorted((float(k), v) for k, v in d.items() if k.replace('.', '', 1).isdigit())
check('F_A: %d certified pairs alpha_k = 0, 0.25, ..., 6.75, all without failure; F_A >= 3.7498 on (0,7]' % len(pairs),
      len(pairs) == 28 and all(v['fails'] == [] for a, v in pairs) and pairs[0][1]['sigma'] == 3.7498
      and all(pairs[i][1]['sigma'] <= pairs[i + 1][1]['sigma'] for i in range(27)))
check('F_A(0+) enclosed (%s); sigma = 3.75 certified NOT a lower bound' % d['FA0'][0], d['neg_375'] is True)
rows = [l.split('#')[0].split() for l in open(os.path.join(C, 'apriori_regions/lamstar.txt'))]
rows = [(float(a), float(b), float(c)) for a, b, c, *_ in [r for r in rows if len(r) >= 3]]
recs = []
for f in ('runs_T0.jsonl', 'runs_T1.jsonl', 'runs_T2.jsonl', 'runs_B.jsonl', 'runs_B2.jsonl'):
    recs += [json.loads(l) for l in open(os.path.join(C, 'apriori_regions/lamstar_runs', f)) if l.strip()]
good = [r for r in recs if r.get('nfails', 1) == 0 and r.get('fails', []) == []]
def covered(lo, hi, L):
    """[lo,hi] is a union of certified runs with Lambda* <= L (runs are on sigma-intervals)"""
    segs = sorted((r['s_lo'], r['s_hi']) for r in good if r['Lstar'] <= L and r['s_hi'] > lo and r['s_lo'] < hi)
    x = lo
    for a, b in segs:
        if a <= x + 1e-12: x = max(x, b)
    return x >= hi - 1e-12
check('lamstar.txt: %d rows, contiguous on [3.7498, 288], no provisional entry' % len(rows),
      rows[0][0] == 3.7498 and rows[-1][1] == 288 and all(rows[i][1] == rows[i + 1][0] for i in range(len(rows) - 1)))
nc = [r for r in rows if not covered(*r)]
check('every row of lamstar.txt is covered by failure-free Arb runs with the same or smaller Lambda* (%d runs)' % len(good), not nc, str(nc[:3]))
code_now = {f: md5(os.path.join(R, 'apriori_regions/lamstar_beta', f)) for f in ('ap_arb.py', 'ap_airy.py', 'ap_airy3.py', 'driver2.py', 'driver3.py', 'driver5.py')}
bad = []
for f in glob.glob(os.path.join(C, 'apriori_regions/apriori_cert/*.json')):
    d = json.load(open(f))
    if any(code_now.get(k) != v for k, v in d['code'].items()): bad.append(os.path.basename(f))
check('145 rectangle certificates: md5 of the code recorded in each = apriori_regions/lamstar_beta/ (validity is checked by sweep.py verify)',
      not bad and len(glob.glob(os.path.join(C, 'apriori_regions/apriori_cert/*.json'))) == 145, ', '.join(bad[:3]))
check('Lambda*(sigma) code: apriori_regions/lamstar_sigma/ap_arb.py has md5 653ed993... (the version reported for these runs)',
      md5(os.path.join(R, 'apriori_regions/lamstar_sigma/ap_arb.py')).startswith('653ed993'))

print('Route 4 (v1.1.0): joint constants of Lemma F5 and Theorem C on [228, 288] (redundant on the logical path, see README)')
R4 = os.path.join(C, 'g1_wall_lemma', 'route4')
if os.path.isdir(R4):
    d = json.load(open(os.path.join(R4, 'M_cert_1_1_0.5.json')))
    nodes = sorted(float(b) for b in d['M_exp'])
    okM = d['all_ok'] is True and d['NJG'] == [1, 1, 0.5] and nodes[0] == 0 and nodes[-1] == 1.0 and len(nodes) == 32
    okM = okM and all(d['M_lower_exp'][b][i] <= d['M_exp'][b][i] for b in d['M_exp'] for i in (0, 1))
    okM = okM and all(d['M_lower_s'][i] <= d['M_s'][i] for i in (0, 1))
    okM = okM and all(r['M'] >= max(r['regions'].values()) for r in d['runs'].values()) and len(d['runs']) == 66
    check('M_u, M_v on 32 nodes beta in [0,1] and M^s_u, M^s_v: enclosures [lower, upper] consistent; M_u(0) <= %.6f, M_v(0) <= %.6f, '
          'M^s = %s' % (d['M_exp']['0'][0], d['M_exp']['0'][1], d['M_s']), okM)
    code_r4 = {f: md5(os.path.join(R, 'g1_wall_lemma', 'route4', f)) for f in d['code']}
    check('md5 of the route-4 code recorded in M_cert = g1_wall_lemma/route4/', code_r4 == d['code'])
    runs = {}
    for f in glob.glob(os.path.join(R4, 'M_runs', 'M_*.jsonl')):
        for l in open(f):
            r = json.loads(l); runs['%s_%s_%r' % (r['pair'], r['weight'], r['beta'])] = r['M']
    check('M_cert values = the 66 recorded runs (M_runs/*.jsonl), rounded up',
          len(runs) == 66 and all(abs(runs[k] - d['runs'][k]['M']) == 0 for k in runs) and
          all({float(x): y for x, y in d['M_exp'].items()}[float(k.split('_')[2])][0 if k[0] == 'u' else 1] >= v for k, v in runs.items() if '_exp_' in k))
    t8 = json.load(open(os.path.join(R4, 'T8_cert_1_1_0.5_s228.json')))
    check('relative-cone table T8 with sigma_min = 228: all %d beta-cells of [0, 0.45] certified' % len(t8['Theta_edge']),
          t8['sigma_min'] == 228.0 and min(a for a, b, th in t8['Theta_edge']) == 0 and max(b for a, b, th in t8['Theta_edge']) == 0.45
          and all(th < 10 for a, b, th in t8['Theta_edge']))
    o = json.load(open(os.path.join(R4, 'out_R4.json')))
    check('Theorem C with route 4 on [228, 288], 0 <= beta <= beta_A: %d boxes, 0 failures, min margin %.3g' % (o['boxes'], o['min_margin']),
          o['fails'] == [] and [float(x) for x in o['NJG']] == [1, 1, 0.5] and o['cA'] == 0.796 and o['coupling'] == 1
          and o['sigma0'] == 228.0 and o['sigma_h'] == 288.0 and o.get('beta_min', 0) == 0)

print('\n%s' % ('ALL CERTIFICATE CHECKS PASSED' if NFAIL == 0 else '%d CHECK(S) FAILED' % NFAIL))
sys.exit(1 if NFAIL else 0)
