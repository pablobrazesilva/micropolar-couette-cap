"""certifica os retangulos de um plano (driver3); em falha de um retangulo com Lambda*, sobe Lambda* em 0.5 e repete.
   uso: python3 run_beta.py plan.json k0 k1 out.jsonl   (retangulos plan[k0:k1])"""
import sys, json, time, os
import driver3 as D3
plan = json.load(open(sys.argv[1])); k0, k1 = int(sys.argv[2]), int(sys.argv[3]); out = sys.argv[4]
os.makedirs('certb', exist_ok=True)
for s_lo, s_hi, b_lo, b_hi, L in plan[k0:k1]:
    L = None if L == 'NONE' else float(L)
    while True:
        tag = 'apriori_s%g-%g_b%g-%g_%s' % (s_lo, s_hi, b_lo, round(b_hi, 4), 'NONE' if L is None else 'L%g' % L)
        t = time.time(); C = D3.Cert(s_lo, s_hi, b_lo, b_hi, L); C.probe = False; C.maxev = 600000
        C.run(1.1)
        res = dict(NJG=[1, 1, 0.5], s_lo=s_lo, s_hi=s_hi, beta_lo=b_lo, beta_hi=b_hi, Lstar=('NONE' if L is None else L),
                   delta_prime='0.05*0.5^(1/3)', P0=C.P0, ratio=1.1, evals=C.nev, boxes_ok=C.nbox, tail_evals=C.tail,
                   skipped=C.skip, fails=[str(f) for f in C.fails[:10]], nfails=len(C.fails), sec=round(time.time() - t, 1),
                   code={f: D3.md5(f) for f in ('ap_arb.py', 'driver2.py', 'driver3.py')}, ok=(len(C.fails) == 0))
        json.dump(res, open('certb/%s.json' % tag, 'w'), indent=1)
        print('RET', json.dumps(dict(tag=tag, ok=res['ok'], evals=res['evals'], sec=res['sec'], fails=res['fails'][:1])), flush=True)
        with open(out, 'a') as f: f.write(json.dumps(dict(res, file='%s.json' % tag)) + '\n')
        if res['ok'] or L is None: break
        L += 0.5
