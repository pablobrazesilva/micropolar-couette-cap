#!/usr/bin/env python3
"""Varredura G2 de K' para (N,J,G) = (1;1;0,5) com o binario `evans` (Arb/FLINT, modelos de Taylor nos parametros).

Subcomandos
  boxes   gera a lista de caixas-base (alpha, sigma) que cobrem K'          -> boxes.txt
  run     processa a fatia --task/--ntasks da lista (subdivisao adaptativa) -> certs/task_XXXX.jsonl
  verify  junta os certificados, confere cada um e a COBERTURA de K'        -> relatorio na saida padrao

K' reduzido (cluster v3) = K' menos as faixas {sigma in [lo,hi], beta >= beta_min} certificadas pelo G1 (g1strip.txt).
K' (relatorio micropolar_G1_arb_v4 / G2 v2, Lema 1.1):  alpha <= alpha0',  F_A(alpha) <= sigma < sigma0'.
  sigma0' = 288 (certificado: crit3 + F8 + rota 3, Arb),  c_A = 0,796 valido para alpha >= 7 (thmA_rig),
  alpha0' = sqrt(288/0,796) = 19,02 -> usamos 19,1.  Cota inferior em sigma usada para DESCARTAR caixas:
  sigma >= c_A alpha^2 para alpha >= 7 (rigoroso);  sigma >= 3,75 = min F_A para todo alpha (ver README, item 'hipoteses').
Lambda*(sigma) vem do ARQUIVO lamstar.txt (entrada do pacote; troque pela tabela certificada).
Cada certificado afirma: para TODO (alpha, sigma) na caixa, D(Lambda) nao tem zeros no retangulo
  [-delta', Lam*] x [-Lam*, S_hi + Lam*]   (delta' = 0,05 g^{1/3}, Lam* = tabela do Teorema 2.1 em sigma_lo),
logo (Teorema 2.1, regiao a priori) nao ha autovalor com Re Lambda >= -delta'  =>  estabilidade espectral na caixa.
"""
import argparse, json, math, os, subprocess, sys, time, glob, hashlib

NJG = ('1', '1', '0.5')
SIGMA_MIN, SIGMA0P, CA, ALPHA_MIN_CA = 3.7498, 288.0, 0.796, 7.0
SIGMA_GRID0 = 3.75          # ancora da grade (compatibilidade com certificados antigos); linha extra [3.7498, 3.75]
ALPHA0P = 19.1
EXE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'evans')

G1STRIP = []          # faixas certificadas pelo G1 (crit3arb): (sigma_lo, sigma_hi, beta_min, arquivo)
def load_g1strip(fn):
    """arquivo g1strip.txt: linhas 'sigma_lo sigma_hi beta_min certificado.json'.  Cada faixa afirma: o criterio do G1
    (crit3 + F8 + rota 3, Arb) vale em {sigma in [lo,hi], beta_min <= beta <= beta_A(sigma)}; acima de beta_A vale o
    Teorema A (alpha >= sqrt(sigma/c_A) >= 7).  Logo nao ha autovalor instavel para sigma in [lo,hi], beta >= beta_min."""
    global G1STRIP
    G1STRIP = []
    if fn and os.path.exists(fn):
        for l in open(fn):
            l = l.split('#')[0].split()
            if len(l) >= 4: G1STRIP.append((float(l[0]), float(l[1]), float(l[2]), l[3]))
    G1STRIP.sort()
def in_g1strip(a1, a2, s1, s2):
    """a caixa inteira esta' na uniao das faixas do G1?  (beta = alpha sigma^{-1/3} >= a1 s^{-1/3}, minimo em s maximo)"""
    if not G1STRIP: return False
    x = s1
    for lo, hi, bmin, fn in G1STRIP:
        if hi <= x or lo > x: continue
        top = min(s2, hi)
        if a1*top**(-1/3) < bmin*(1 + 1e-12): return False
        x = hi
        if x >= s2: return True
    return False

FAPAIRS = []      # pares certificados (alpha_k, sigma_k): F_A >= sigma_k em [alpha_k, 7]  (fa_cert.json, G2 apriori Arb)
def load_fa(fn='fa_cert.json'):
    global FAPAIRS
    fn = fn if os.path.isabs(fn) else os.path.join(os.path.dirname(os.path.abspath(__file__)), fn)
    d = json.load(open(fn)); FAPAIRS = []
    for k, v in d.items():
        try: a = float(k)
        except ValueError: continue
        if isinstance(v, dict) and v.get('fails') == [] and 'sigma' in v: FAPAIRS.append((a, float(v['sigma'])))
    FAPAIRS.sort()
def fa_lower(a1, a2):
    """cota inferior rigorosa de F_A em [a1, a2]"""
    if a1 >= ALPHA_MIN_CA: return CA*a1*a1
    b = max([sk for ak, sk in FAPAIRS if ak <= a1] + [0.0])
    if a2 > ALPHA_MIN_CA: b = min(b, CA*ALPHA_MIN_CA**2)
    return b

def needed(a1, a2, s1, s2):
    """a caixa intersecta K' reduzido?  (descarta so' com argumentos rigorosos)"""
    if s1 >= SIGMA0P or s2 < SIGMA_MIN or a1 > ALPHA0P: return False
    if not FAPAIRS: load_fa()
    if s2 < fa_lower(a1, a2): return False                          # Teorema A / pares F_A: sigma < F_A(alpha)
    if lam_beta_none(a1, a2, s1, s2): return False                  # regiao a priori vazia (lamstar_beta.txt, NONE)
    if in_g1strip(a1, a2, s1, s2): return False                      # G1 em faixas de beta (cluster v3)
    return True

def check_g1strip():
    """confere cada faixa contra o certificado Arb do driver G1 (out_*.json)"""
    bad = []
    for lo, hi, bmin, fn in G1STRIP:
        try:
            d = json.load(open(fn if os.path.isabs(fn) else os.path.join(os.path.dirname(os.path.abspath(__file__)), fn)))
            ok = (d['fails'] == [] and [float(x) for x in d['NJG']] == [1.0, 1.0, 0.5] and d['sigma0'] <= lo + 1e-9
                  and d['sigma_h'] >= hi - 1e-9 and d.get('beta_min', 1e9) <= bmin + 1e-12 and abs(d['cA'] - CA) < 1e-12
                  and d.get('coupling', 1) == 1 and lo >= CA*ALPHA_MIN_CA**2)
        except Exception as e:
            ok = False
        if not ok: bad.append((lo, hi, bmin, fn))
    return bad

def base_boxes(da=0.5, dls=0.05):
    out = []
    ns = int(math.ceil(math.log(SIGMA0P/SIGMA_GRID0)/dls)); na = int(math.ceil(ALPHA0P/da))
    rows = [(SIGMA_MIN, SIGMA_GRID0)] + [(SIGMA_GRID0*math.exp(i*dls), min(SIGMA0P, SIGMA_GRID0*math.exp((i + 1)*dls))) for i in range(ns)]
    for (s1, s2) in rows:
        for j in range(na):
            a1, a2 = j*da, min(ALPHA0P, (j + 1)*da)
            if needed(a1, a2, s1, s2): out.append((a1, a2, s1, s2))
    return out

LAMSTAR = None; LAMSTAR_SHA = None
def load_lamstar(fn):
    """arquivo de entrada Lambda*(sigma): linhas 'sigma_lo sigma_hi Lstar [flag]'"""
    global LAMSTAR, LAMSTAR_SHA
    rows = []
    for l in open(fn):
        l = l.split('#')[0].split()
        if len(l) >= 3: rows.append((float(l[0]), float(l[1]), float(l[2]), l[3] if len(l) > 3 else ''))
    LAMSTAR = rows; LAMSTAR_SHA = hashlib.sha256(open(fn, 'rb').read()).hexdigest()[:16]
def lam_star(s1, s2):
    """max Lambda* sobre as faixas que intersectam [s1, s2]; None se [s1,s2] nao estiver coberto"""
    v = None; cov = []
    for lo, hi, L, fl in LAMSTAR:
        if hi > s1 and lo < s2: v = L if v is None else max(v, L); cov.append((lo, hi))
    cov.sort(); x = s1
    for lo, hi in cov:
        if lo > x + 1e-12: return None
        x = max(x, hi)
    return v if x >= s2 - 1e-12 else None
def lam_flags(s1, s2): return sorted(set(fl for lo, hi, L, fl in LAMSTAR if hi > s1 and lo < s2 and fl))

# ---- regiao a priori dependente de beta (formato combinado; arquivo do agente matematico) ----
LAMB = None; LAMB_SHA = None; LAMB_FILE = None
def load_lamstar_beta(fn):
    """lamstar_beta.txt: 'sigma_lo sigma_hi beta_lo beta_hi Lamstar cert.json', Lamstar numero ou NONE; '#' comenta.
    Retangulos FECHADOS; sobreposicao permitida (vale o menor; NONE < qualquer numero).  Ausente -> None (fallback)."""
    global LAMB, LAMB_SHA, LAMB_FILE
    LAMB = None; LAMB_SHA = None; LAMB_FILE = None
    if not fn or not os.path.exists(fn): return
    rows = []
    for l in open(fn):
        t = l.split('#')[0].split()
        if len(t) < 6: continue
        L = None if t[4].upper() == 'NONE' else float(t[4])
        rows.append((float(t[0]), float(t[1]), float(t[2]), float(t[3]), L, t[5]))
    LAMB = rows; LAMB_FILE = fn; LAMB_SHA = hashlib.sha256(open(fn, 'rb').read()).hexdigest()[:16]
def _lamb_box(a1, a2, s1, s2):
    """(todo_none, Lmax, coberto) para a imagem da caixa em (sigma, beta), beta = alpha sigma^{-1/3}.
    Fatias em sigma nas fronteiras dos retangulos; em cada fatia [u,v], beta in [a1 v^{-1/3}, a2 u^{-1/3}] (superconjunto).
    Em cada subintervalo elementar de beta, o valor e' o minimo dos retangulos que o contem (NONE = sem autovalores);
    nos pontos de fronteira o minimo so' pode ser menor (retangulos fechados), entao o maximo sobre os subintervalos basta."""
    cuts = sorted(set([s1, s2] + [x for r in LAMB for x in (r[0], r[1]) if s1 < x < s2]))
    allnone = True; Lmax = None; covered = True; anynone = False
    for u, v in zip(cuts[:-1], cuts[1:]):
        b1, b2 = a1*v**(-1/3)*(1 - 1e-12), a2*u**(-1/3)*(1 + 1e-12)
        rs = [r for r in LAMB if r[0] <= u and r[1] >= v and r[3] >= b1 and r[2] <= b2]
        bc = sorted(set([b1, b2] + [x for r in rs for x in (r[2], r[3]) if b1 < x < b2]))
        for x, y in zip(bc[:-1], bc[1:]):
            cov = [r for r in rs if r[2] <= x and r[3] >= y]
            if not cov: covered = False; allnone = False; continue
            if any(r[4] is None for r in cov): anynone = True; continue
            allnone = False
            L = min(r[4] for r in cov); Lmax = L if Lmax is None else max(Lmax, L)
    _lamb_box.anynone = anynone
    return allnone and covered, Lmax, covered
def lam_beta_partial(a1, a2, s1, s2):
    """a caixa toca algum retangulo NONE sem estar inteira em NONE?"""
    if LAMB is None: return False
    allnone = _lamb_box(a1, a2, s1, s2)[0]
    return (not allnone) and _lamb_box.anynone
def lam_beta_none(a1, a2, s1, s2):
    if LAMB is None: return False
    return _lamb_box(a1, a2, s1, s2)[0]
def lam_box(a1, a2, s1, s2):
    """Lambda* para a caixa: tabela por beta onde cobre; nas partes nao cobertas, a tabela so' em sigma (valida em K')."""
    Ls = lam_star(s1, s2) if LAMSTAR is not None else None
    if LAMB is None: return Ls
    allnone, Lb, cov = _lamb_box(a1, a2, s1, s2)
    if allnone: return 0.0
    if cov: return Lb if Lb is not None else 0.0
    if Ls is None: return None
    return max(Ls, Lb) if Lb is not None else Ls
def check_lamb_certs():
    bad = []
    if LAMB is None: return bad
    base = os.path.dirname(os.path.abspath(LAMB_FILE))
    for r in LAMB:
        try:
            d = json.load(open(os.path.join(base, r[5])))
            Lc = d.get('Lstar'); Lc = None if (isinstance(Lc, str) and Lc.upper() == 'NONE') else float(Lc)
            ok = (d.get('fails') == [] and d.get('ok') is True and [float(x) for x in d['NJG']] == [1.0, 1.0, 0.5]
                  and float(d['s_lo']) <= r[0] + 1e-12 and float(d['s_hi']) >= r[1] - 1e-12
                  and float(d['beta_lo']) <= r[2] + 1e-12 and float(d['beta_hi']) >= r[3] - 1e-12
                  and ((r[4] is None and Lc is None) or (r[4] is not None and Lc is not None and Lc <= r[4] + 1e-12)))
        except Exception:
            ok = False
        if not ok: bad.append(r)
    return bad

RUN_MODE = 'contorno'
def run_box(a1, a2, s1, s2, timeout):
    Ls = lam_box(a1, a2, s1, s2)
    if Ls is None: return dict(ok=0, error='Lambda* nao coberto pelo arquivo', alpha=[a1, a2], sigma=[s1, s2])
    cmd = [EXE, 'box', *NJG, repr(a1), repr(a2), repr(s1), repr(s2), '6', '0.25', '16', repr(Ls)]
    env = dict(os.environ)
    if RUN_MODE == 'continuacao':
        r = row_of((a1, a2, s1, s2))
        if r is not None and row_anchor(r): env['EVANS_MODE'] = 'line'
    t = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
        d = json.loads(p.stdout.strip().splitlines()[-1])
    except Exception as e:
        d = dict(ok=0, error=str(e)[:200], alpha=[a1, a2], sigma=[s1, s2])
    d['wall'] = time.time() - t; d['lamstar_sha'] = LAMSTAR_SHA; d['lamstar_beta_sha'] = LAMB_SHA; d['lamstar_flags'] = lam_flags(s1, s2)
    return d

def children(a1, a2, s1, s2):
    """subdivisao 2x2 canonica (alpha linear, sigma geometrico) -- a MESMA usada em run e verify"""
    am, sm = 0.5*(a1 + a2), math.sqrt(s1*s2)
    return [(a1, am, s1, sm), (am, a2, s1, sm), (a1, am, sm, s2), (am, a2, sm, s2)]
def band_of(b): return 'A' if b[2] < 20 else ('B' if b[2] < 60 else 'C')
def beta_range(b): return b[0]*b[3]**(-1/3), b[1]*b[2]**(-1/3)
DEPTH_NONE, DEPTH_CMID, CMID = 3, 2, (0.5, 1.5)
def grid_rows(dls=0.05):
    ns = int(math.ceil(math.log(SIGMA0P/SIGMA_GRID0)/dls))
    return [(SIGMA_MIN, SIGMA_GRID0)] + [(SIGMA_GRID0*math.exp(i*dls), min(SIGMA0P, SIGMA_GRID0*math.exp((i + 1)*dls))) for i in range(ns)]
def row_of(b, dls=0.05):
    for r in grid_rows(dls):
        if r[0] <= b[2] + 1e-12 and b[3] <= r[1] + 1e-12: return r
    return None
def row_anchor(row, da=0.5, cert=None):
    """modo continuacao (Teorema 3.2 do relatorio apriori-beta v2, na forma da secao 3 do cluster v5): na linha de grade
    [S1,S2], percorre as caixas-base em alpha crescente a partir de 0; devolve True se chega a uma caixa INTEIRA em NONE
    antes de sair do K' reduzido ou da cobertura dos retangulos a priori, e (se cert for dado) se toda caixa necessaria
    no caminho estiver coberta por certificados (contorno ou reta)."""
    if LAMB is None: return False
    j = 0
    while j*da <= ALPHA0P + da:
        b = (j*da, min(ALPHA0P, (j + 1)*da), row[0], row[1]); j += 1
        if lam_beta_none(*b): return True
        if not needed(*b): return False
        if not _lamb_box(*b)[2]: return False
        if cert is not None and not covered_by(cert, b): return False
    return False

def box_cost(b):
    """custo relativo de uma caixa ~ comprimento do contorno 2(S_hi + 2 L*) + 2 L* (medido: tempo ~ no. de segmentos)"""
    L = lam_box(*b) or 0.0
    return 2*(2*b[3]**(1/3) + 2*L) + 2*L + 0.1
def leaves(b, level=0):
    """pre-subdivisao (folhas da arvore 2x2 canonica):
    (i) caixa que cruza a fronteira NONE: divide-se (ate DEPTH_NONE niveis) SO' se a soma dos custos das filhas que
        sobram (as inteiras em NONE saem por needed) for menor que o custo da caixa -- o custo cai com o perimetro do
        contorno (Lambda* menor) e com as filhas que saem;
    (ii) faixa C com beta em [0,5; 1,5]: divide-se ate DEPTH_CMID (resto do modelo de Taylor, medido na v2-v4)."""
    if not needed(*b): return [], 0.0
    bl, bh = beta_range(b)
    if level < DEPTH_CMID and band_of(b) == 'C' and bh >= CMID[0] and bl <= CMID[1]:
        out, cost = [], 0.0
        for c in children(*b):
            o, k = leaves(c, level + 1); out += o; cost += k
        return out, cost
    me = box_cost(b)
    if level < 1 and band_of(b) in 'BC' and lam_beta_partial(*b):
        # medido (v5): caixas-base de B/C que cruzam a fronteira NONE falham na 1a tentativa (50-200 s perdidos) e
        # passam no nivel 1 -> divide-se sempre um nivel; as filhas inteiras em NONE saem
        out, cost = [], 0.0
        for c in children(*b):
            o, k = leaves(c, level + 1); out += o; cost += k
        return out, cost
    if level < DEPTH_NONE and lam_beta_partial(*b):
        out, cost = [], 0.0
        for c in children(*b):
            o, k = leaves(c, level + 1); out += o; cost += k
        if cost < me: return out, cost
    return [(b, level)], me

def cmd_boxes(args):
    load_g1strip(args.g1strip); load_lamstar_beta(args.lamstar_beta); load_fa()
    bx = base_boxes(args.da, args.dls)
    with open(args.out, 'w') as f:
        for b in bx: f.write('%r %r %r %r\n' % b)
    print('%d caixas-base -> %s' % (len(bx), args.out))
    # lista de TAREFAS (folhas da pre-subdivisao), por faixa
    tk = {'A': [], 'B': [], 'C': []}
    for b in bx:
        for (c, lev) in leaves(b)[0]: tk[band_of(b)].append(c)
    stem = args.out[:-4] if args.out.endswith('.txt') else args.out
    for k, v in tk.items():
        with open('tasks_%s.txt' % k if stem == 'boxes' else stem + '_tasks_%s.txt' % k, 'w') as f:
            for c in v: f.write('%r %r %r %r\n' % c)
    print('tarefas (folhas pre-subdivididas): A %d, B %d, C %d -> tasks_A/B/C.txt' % (len(tk['A']), len(tk['B']), len(tk['C'])))

def cmd_split(args):
    bx = [tuple(map(float, l.split())) for l in open(args.boxes) if l.strip()]
    bands = (('A', 0, 20), ('B', 20, 60), ('C', 60, 1e9))
    for name, lo, hi in bands:
        sel = [b for b in bx if lo <= b[2] < hi]
        with open(args.boxes.replace('.txt', '_%s.txt' % name), 'w') as f:
            for b in sel: f.write('%r %r %r %r\n' % b)
        print('faixa %s (%g <= sigma_lo < %g): %d caixas-base' % (name, lo, hi, len(sel)))

def cert_valid(d):
    """certificado ainda valido com as entradas ATUAIS (Lambda* usado >= exigido agora; terna; ok)"""
    key = (d['alpha'][0], d['alpha'][1], d['sigma'][0], d['sigma'][1])
    Lreq = lam_box(*key)
    base = (d.get('ok') == 1 and list(d.get('NJG', [])) == list(NJG) and Lreq is not None and d.get('Lstar', 0) >= Lreq - 1e-12)
    if d.get('mode') == 'line': return base and d.get('line_ok') == 1, key
    return base and d.get('count_certified') == 1, key
def valid_certs(files, kind='all'):
    """kind: 'all' (contorno + reta), 'contour', 'line'"""
    out = {}
    for fn in files:
        if not os.path.exists(fn): continue
        for l in open(fn):
            try: d = json.loads(l)
            except Exception: continue
            ok, key = cert_valid(d)
            m = d.get('mode', 'contour')
            if ok and (kind == 'all' or kind == m): out[key] = d
    return out
# ---------------- modo continuacao: ancoragem por componentes conexas (Prop. 5.13 do artigo, condicoes V1/V2) ----------------
TOL_A, TOL_S = 1e-12, 1e-12
def _touch(b, c):
    """adjacencia de duas caixas fechadas em (alpha, sigma): interiores que se intersectam, OU face comum com interseccao
    de COMPRIMENTO POSITIVO.  Contato so' por canto NAO conta (escolha conservadora: evita degenerescencias de ponto
    flutuante; todo caminho por aresta e' admissivel na Prop. 5.13)."""
    ao = min(b[1], c[1]) - max(b[0], c[0]); so = min(b[3], c[3]) - max(b[2], c[2])
    ta = TOL_A*max(1.0, abs(b[1]), abs(c[1])); ts = TOL_S*max(b[3], c[3])
    if ao > ta and so > ts: return True                      # sobreposicao de area positiva
    if abs(ao) <= ta and so > ts: return True                # face vertical comum (alpha = const)
    if abs(so) <= ts and ao > ta: return True                # face horizontal comum (sigma = const)
    return False
def anchor_components(nodes):
    """nodes: lista de (a1, a2, s1, s2, kind), kind in {'line', 'contour', 'none'}.  Union-find sobre a adjacencia _touch.
    Devolve o conjunto de indices de nos 'line' cuja componente contem uma ancora ('contour' ou 'none').
    Caixas descartadas so' pelo Teorema A (ou pelo G1) NAO sao nos: nao servem de ponte."""
    n = len(nodes); par = list(range(n))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    order = sorted(range(n), key=lambda i: nodes[i][2])          # varredura em sigma
    active = []
    for i in order:
        b = nodes[i]
        active = [j for j in active if nodes[j][3] >= b[2] - TOL_S*b[2]]
        for j in active:
            if _touch(b, nodes[j]):
                ri, rj = find(i), find(j)
                if ri != rj: par[ri] = rj
        active.append(i)
    anch = set(find(i) for i in range(n) if nodes[i][4] in ('contour', 'none'))
    return set(i for i in range(n) if nodes[i][4] == 'line' and find(i) in anch)
def none_nodes(depth=4, da=0.5, dls=0.05):
    """sub-regioes INTEIRAMENTE em retangulos NONE (pior canto), como folhas da arvore 2x2 canonica sobre a grade toda"""
    out = []
    if LAMB is None: return out
    def rec(b, d):
        if lam_beta_none(*b): out.append(b + ('none',)); return
        if d >= depth or not lam_beta_partial(*b): return
        for c in children(*b): rec(c, d + 1)
    na = int(math.ceil((ALPHA0P + da)/da))
    for (s1, s2) in grid_rows(dls):
        for j in range(na): rec((j*da, (j + 1)*da, s1, s2), 0)
    return out
def effective_certs(files):
    """certificados que valem: os de contorno, mais os de reta que satisfazem V1 (caixa inteira coberta pelos retangulos
    certificados de lamstar_beta.txt, Lambda* usado >= exigido) e V2 (mesma componente conexa, por adjacencia de face,
    que uma caixa NONE ou de contorno certificada)."""
    full = valid_certs(files, 'contour'); line = valid_certs(files, 'line')
    line_v1 = {k: d for k, d in line.items() if LAMB is not None and _lamb_box(*k)[2]}
    nodes = [k + ('contour',) for k in full] + [k + ('line',) for k in line_v1] + none_nodes()
    ok = anchor_components(nodes)
    acc = set(nodes[i][:4] for i in ok)
    eff = dict(full)
    for k, d in line_v1.items():
        if k in acc: eff[k] = d
    return eff, full, line, acc

def covered_by(cert, b, depth=0):
    """b coberta: fora do K' reduzido, certificada, ou suas 4 filhas cobertas (recursivo)"""
    if not needed(*b): return True
    if b in cert: return True
    if depth >= 12: return False
    if not any(k for k in cert if k[0] >= b[0] and k[1] <= b[1] and k[2] >= b[2] and k[3] <= b[3]): return False
    return all(covered_by(cert, c, depth + 1) for c in children(*b))

def cmd_todo(args):
    """lista de tarefas ainda NAO cobertas por certificados validos (apos mudar lamstar_beta.txt, etc.)"""
    load_lamstar(args.lamstar); load_g1strip(args.g1strip); load_lamstar_beta(args.lamstar_beta); load_fa()
    cert = effective_certs(glob.glob(args.certs))[0]
    tot = 0; left = []
    for fn in args.tasks:
        for l in open(fn):
            if not l.strip(): continue
            b = tuple(map(float, l.split())); tot += 1
            if not covered_by(cert, b): left.append(b)
    with open(args.out, 'w') as f:
        for b in left: f.write('%r %r %r %r\n' % b)
    print('certificados validos: %d; tarefas: %d; ainda a fazer: %d -> %s' % (len(cert), tot, len(left), args.out))

def cmd_run(args):
    global RUN_MODE; RUN_MODE = args.modo
    load_lamstar(args.lamstar); load_g1strip(args.g1strip); load_lamstar_beta(args.lamstar_beta); load_fa()
    bx = [tuple(map(float, l.split())) for l in open(args.boxes) if l.strip()]
    mine = bx[args.task::args.ntasks]
    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    done = valid_certs([args.out] + glob.glob(args.certs))      # retomada idempotente: so' certificados AINDA validos
    with open(args.out, 'a') as f:
        for b in mine:
            stack = [(b, 0)]
            while stack:
                (a1, a2, s1, s2), depth = stack.pop()
                if not needed(a1, a2, s1, s2) or covered_by(done, (a1, a2, s1, s2)): continue
                d = run_box(a1, a2, s1, s2, args.timeout); d['depth'] = depth; d['base'] = list(b)
                f.write(json.dumps(d) + '\n'); f.flush()
                if d.get('ok') and (d.get('count_certified') or d.get('line_ok')):
                    if d.get('count', 0) != 0: print('ATENCAO: %d zero(s) na caixa %s' % (d['count'], (a1, a2, s1, s2)), flush=True)
                    continue
                if depth >= args.maxdepth: print('FALHA (profundidade maxima) %s' % ((a1, a2, s1, s2),), flush=True); continue
                stack += [(c, depth + 1) for c in children(a1, a2, s1, s2)]

def cmd_verify(args):
    load_lamstar(args.lamstar); load_g1strip(args.g1strip); load_lamstar_beta(args.lamstar_beta); load_fa()
    badl = check_lamb_certs()
    print('regiao a priori por beta: %s (%s retangulos, sha256 %s), certificados invalidos: %d'
          % (LAMB_FILE or 'ausente -> fallback lamstar.txt', len(LAMB) if LAMB else 0, LAMB_SHA, len(badl)))
    for b in badl[:10]: print('  RETANGULO SEM CERTIFICADO VALIDO:', b)
    badstrip = check_g1strip()
    print('faixas G1 (g1strip): %d, invalidas: %d' % (len(G1STRIP), len(badstrip)))
    for b in badstrip: print('  FAIXA G1 INVALIDA:', b)
    bx = [tuple(map(float, l.split())) for l in open(args.boxes) if l.strip()]
    cert = {}; nbad = 0; zeros = []; tot_t = 0.0; nlines = 0; provisional = set()
    nobsolete = 0
    for fn in glob.glob(args.certs):
        for l in open(fn):
            d = json.loads(l); nlines += 1; tot_t += d.get('seconds', d.get('wall', 0))
            good, key = cert_valid(d)
            if not good and d.get('ok') == 1 and d.get('count_certified') == 1: nobsolete += 1
            if good and lam_flags(key[2], key[3]): provisional.add(key)
            if good:
                cert[key] = d
                if d.get('count', 0) != 0: zeros.append(key)
    # cobertura: cada caixa-base e' coberta se ela ou (recursivamente) seus 4 filhos estao certificados ou fora de K'
    eff, full, line, okrows = effective_certs(glob.glob(args.certs))
    print('certificados de contorno: %d; de reta (continuacao): %d, aceitos (V1+V2, componente com ancora): %d'
          % (len(full), len(line), sum(1 for k in line if k in eff)))
    cert = eff
    zeros = [k for k in zeros if k in cert]
    missing = [b for b in bx if not covered_by(cert, b)]
    print('certificados obsoletos (Lambda* insuficiente para as entradas atuais), ignorados: %d' % nobsolete)
    # a lista-base cobre K'?  (faixas contiguas em sigma e alpha)
    print('linhas lidas: %d   certificados validos: %d   tempo total (CPU dos certificados): %.1f h' % (nlines, len(cert), tot_t/3600))
    print('caixas-base: %d   descobertas: %d   caixas com zeros (instabilidade!): %d' % (len(bx), len(missing), len(zeros)))
    for m in missing[:20]: print('  descoberta:', m)
    for z in zeros[:20]: print('  ZERO em:', z)
    ref = base_boxes(args.da, args.dls)
    same = sorted(ref) == sorted(bx)
    print('lista-base identica a gerada por base_boxes(%g, %g): %s' % (args.da, args.dls, same))
    print('arquivo Lambda*: %s (sha256 %s); certificados que usam faixas PROVISORIAS: %d' % (args.lamstar, LAMSTAR_SHA, len(provisional)))
    ok = not missing and not zeros and same and not badstrip and not badl
    print('RESULTADO:', ('K\' COBERTO, 0 autovalores instaveis' + (' (CONDICIONAL: ha faixas PROVISORIAS de Lambda*)' if provisional else ''))
          if ok else 'INCOMPLETO')
    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    b = sp.add_parser('boxes'); b.add_argument('--lamstar-beta', dest='lamstar_beta', default='lamstar_beta.txt'); b.add_argument('--g1strip', default='g1strip.txt'); b.add_argument('--da', type=float, default=0.5); b.add_argument('--dls', type=float, default=0.05); b.add_argument('--out', default='boxes.txt')
    r = sp.add_parser('run'); r.add_argument('--modo', choices=['contorno', 'continuacao'], default='contorno'); r.add_argument('--lamstar-beta', dest='lamstar_beta', default='lamstar_beta.txt'); r.add_argument('--certs', default='certs/*.jsonl'); r.add_argument('--boxes', default='boxes.txt'); r.add_argument('--task', type=int, required=True); r.add_argument('--ntasks', type=int, required=True)
    r.add_argument('--lamstar', default='lamstar.txt'); r.add_argument('--g1strip', default='g1strip.txt'); r.add_argument('--out', required=True); r.add_argument('--timeout', type=float, default=7200); r.add_argument('--maxdepth', type=int, default=4)
    v = sp.add_parser('verify'); v.add_argument('--lamstar-beta', dest='lamstar_beta', default='lamstar_beta.txt'); v.add_argument('--boxes', default='boxes.txt'); v.add_argument('--certs', default='certs/*.jsonl')
    v.add_argument('--lamstar', default='lamstar.txt'); v.add_argument('--g1strip', default='g1strip.txt'); v.add_argument('--da', type=float, default=0.5); v.add_argument('--dls', type=float, default=0.05)
    sp2 = sp.add_parser('split'); sp2.add_argument('--boxes', default='boxes.txt')
    t = sp.add_parser('todo'); t.add_argument('--lamstar', default='lamstar.txt'); t.add_argument('--lamstar-beta', dest='lamstar_beta', default='lamstar_beta.txt')
    t.add_argument('--g1strip', default='g1strip.txt'); t.add_argument('--certs', default='certs/*.jsonl'); t.add_argument('--out', default='todo.txt')
    t.add_argument('tasks', nargs='+')
    a = ap.parse_args()
    {'boxes': cmd_boxes, 'run': cmd_run, 'verify': cmd_verify, 'split': cmd_split, 'todo': cmd_todo}[a.cmd](a)
