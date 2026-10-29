"""Teste da idempotencia em relacao a lamstar_beta.txt: com os certificados finais (originalmente: do piloto),
(1) entradas atuais -> 'todo' exclui as tarefas ja' cobertas;
(2) uma copia de lamstar_beta.txt com TODOS os Lambda* numericos dobrados -> os certificados com Lambda* menor viram
    obsoletos e as tarefas correspondentes voltam para 'todo'.  (Arquivo de teste temporario; nao altera o original.)"""
import os, sys, subprocess, tempfile, glob
here = os.path.dirname(os.path.abspath(__file__)); pkg = os.path.join(here, '..')
certs = os.path.join(pkg, '..', 'certificates', 'g2_evans', 'certs_final', '*.jsonl')   # [repo] the pilot certificates are not distributed; use the final ones
tasks = [os.path.join(pkg, 'tasks_%s.txt' % k) for k in 'ABC']
def todo(lb):
    out = tempfile.mktemp(suffix='.txt')
    r = subprocess.run([sys.executable, os.path.join(pkg, 'sweep.py'), 'todo', '--lamstar', os.path.join(pkg, 'lamstar.txt'),
                        '--lamstar-beta', lb, '--g1strip', os.path.join(pkg, 'g1strip.txt'), '--certs', certs, '--out', out] + tasks,
                       capture_output=True, text=True, cwd=pkg)
    print(r.stdout.strip()); n = sum(1 for l in open(out) if l.strip()); os.remove(out); return n
n0 = todo(os.path.join(pkg, 'lamstar_beta.txt'))
tmp = os.path.join(pkg, 'lamstar_beta_TESTE_dobrado.txt')
with open(os.path.join(pkg, 'lamstar_beta.txt')) as f, open(tmp, 'w') as g:
    for l in f:
        t = l.split('#')[0].split()
        if len(t) >= 6 and t[4].upper() != 'NONE': t[4] = repr(2*float(t[4])); l = ' '.join(t) + '\n'
        g.write(l)
n1 = todo(tmp); os.remove(tmp)
ntot = sum(1 for fn in tasks for l in open(fn) if l.strip())
print('tarefas %d; a fazer com entradas atuais %d; com Lambda* dobrado %d' % (ntot, n0, n1))
ok = n0 < ntot and n1 > n0
print('OK' if ok else 'FALHA'); sys.exit(0 if ok else 1)
