"""roda driver5 (NONE por continuacao) para uma lista de retangulos [s_lo, s_hi, b_lo, b_hi]; em falha sobe b_lo em 0,05."""
import sys, json, time, os
import driver5 as D5
rects = json.load(open(sys.argv[1])); out = sys.argv[2]
os.makedirs('certc', exist_ok=True)
for s_lo, s_hi, b_lo, b_hi in rects:
    while b_lo < b_hi - 1e-9:
        tag = 'apriori_cont_s%g-%g_b%g-%g_NONE' % (s_lo, s_hi, b_lo, b_hi)
        t = time.time(); C = D5.LineCert(s_lo, s_hi, b_lo, b_hi); C.run(verbose=False)
        res = dict(NJG=[1, 1, 0.5], s_lo=s_lo, s_hi=s_hi, beta_lo=b_lo, beta_hi=b_hi, Lstar='NONE',
                   metodo='continuacao: reta Re Lambda = -delta\' + cauda, ancorada no retangulo NONE adjacente em beta = beta_hi',
                   delta_prime='0.05*0.5^(1/3)', P0=C.P0, evals=C.nev, LG_ok=C.nLG, airy_evals=C.nAE, tail_evals=C.tail,
                   fails=[str(f) for f in C.fails[:10]], nfails=len(C.fails), sec=round(time.time() - t, 1),
                   code={f: D5.md5(f) for f in ('ap_arb.py', 'ap_airy.py', 'ap_airy3.py', 'driver2.py', 'driver5.py')},
                   ok=(len(C.fails) == 0))
        json.dump(res, open('certc/%s.json' % tag, 'w'), indent=1)
        print('RET', json.dumps(dict(tag=tag, ok=res['ok'], evals=res['evals'], airy=res['airy_evals'], sec=res['sec'], fails=res['fails'][:1])), flush=True)
        with open(out, 'a') as f: f.write(json.dumps(dict(res, file=tag + '.json')) + '\n')
        if res['ok']: break
        b_lo = round(b_lo + 0.05, 2)
