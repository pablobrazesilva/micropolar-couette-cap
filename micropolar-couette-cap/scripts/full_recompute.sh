#!/bin/sh
# Level 3: full recomputation of every computational link of the logical path, component by component.
#   sh scripts/full_recompute.sh [g0] [consts] [pi] [g1] [thmA] [lamstar] [lambeta] [route4] [g2]  (default: all but route4, g2)
# Requires FLINT 3.0.1, python-flint 0.9.0, numpy, scipy, mpmath.  New outputs go to $CAP_WORK (default ./work);
# certificates/ is never modified.  Compare the outputs with certificates/ (same file names), or rerun
# scripts/check_certificates.py on a copy of the repository in which certificates/ has been replaced by the new files.
# CPU times of the original runs (2 cores): g0 1 min, consts 30 min, pi seconds, g1 30 min, thmA seconds,
# lamstar 1 h, lambeta 4 h, g2 28.5 h.
set -e
PY0="${PYTHON:-python3}"; case "$PY0" in */*) PY0="$(cd "$(dirname "$PY0")" && pwd)/$(basename "$PY0")";; esac   # absolute path, resolved before cd
cd "$(dirname "$0")/.."
R="$(pwd)"; PY="$PY0"
export CAP_WORK="${CAP_WORK:-$R/work}"; mkdir -p "$CAP_WORK"
STEPS="${*:-g0 consts pi g1 thmA lamstar lambeta}"
wcopy() { mkdir -p "$CAP_WORK/$1"; cp -R "$R/$1/." "$CAP_WORK/$1/"; cd "$CAP_WORK/$1"; }   # a priori scripts run in a copy

for s in $STEPS; do case $s in
g0)  # Appendix A: cone of G at eta = 0.05 and 0.1 (all rho >= 0), then the negative controls
  $PY g0_sector_lemma/run_W.py 0.05 -0.0836 0.54937 0.6083 20 140 40 C1_005
  $PY g0_sector_lemma/run_W.py 0.1  -0.0837 0.56137 0.5827 20 140 40 C1_010
  $PY g0_sector_lemma/audit/negctl.py ;;
consts)  # Appendices B and D (eta = 0.05).  The later steps read the RECORDED line_0.05.json, tables and T8.
  $PY g0_constants/certline.py 0.05 20 30 0.00390625
  $PY g0_constants/refine.py 0.05 u 0.44452699;  $PY g0_constants/refine.py 0.05 ut 0.35346270
  $PY g0_constants/refine.py 0.05 ut 0.35346272 || echo "negative control: not certified (expected)"
  $PY g0_constants/pbound.py 0.05;  $PY g0_constants/ptable.py
  $PY g0_constants/tp_cert.py;      $PY g0_constants/cert_T36.py
  $PY g0_constants/cert_J2.py 1 1 0.5
  $PY g0_constants/cert_T8.py 1 1 0.5 240 3.3;  $PY g0_constants/t8_patch.py 1 1 0.5 0.04 0.3 1.3
  $PY g0_constants/t8_sub.py 1 1 0.5 215 0.45 3
  $PY g0_constants/t8_sub.py 1 1 0.5 215 1.7 3 0.25   # second pass on the two coarse cells (its tolerance was not logged)
  $PY g0_constants/controls.py ;;                      # cone negative controls and point lower bounds of P (needs scipy)
pi)
  $PY g1_wall_lemma/pi_tails/cert_pi_arb.py
  $PY g1_wall_lemma/pi_tails/controles.py;  $PY g1_wall_lemma/pi_tails/audit_T36_caudas.py ;;
g1)  # Appendix D: Theorem C criterion, sigma >= 288, all 0 <= beta <= beta_A; tail sigma >= 3e5; bands below 288
  G1_CRIT=crit3arb G1_MAXEVAL=2000 $PY g1_wall_lemma/driver.py 1 1 0.5 0.796 288 300 C4R2b 1
  G1_CRIT=crit3arb G1_MAXEVAL=3000 $PY g1_wall_lemma/driver.py 1 1 0.5 0.796 300 650 C4R2 1
  G1_CRIT=crit3arb G1_MAXEVAL=20000 $PY g1_wall_lemma/driver.py 1 1 0.5 0.796 650 1100 C3R2 1
  G1_CRIT=crit3arb G1_ROUTE3=0 $PY g1_wall_lemma/driver.py 1 1 0.5 0.796 1100 1900 C3A2 1
  $PY g1_wall_lemma/driver.py 1 1 0.5 0.796 1900 3e5 A2 1
  $PY g1_wall_lemma/tail.py 1 1 0.5 0.796 3e5
  for b in "280 288 0.05" "265 280 0.1" "250 265 0.2" "240 250 0.4" "232 240 0.45" "225 232 0.5" "215 225 0.6"; do
    set -- $b; G1_CRIT=crit3arb G1_MAXEVAL=3000 G1_BMIN=$3 $PY g1_wall_lemma/driver.py 1 1 0.5 0.796 $1 $2 S_$1 1; done
  $PY g1_wall_lemma/negctl_crit3.py 1 1 0.5 260;  $PY g1_wall_lemma/negctl_crit3.py 1 1 0.5 275
  $PY g1_wall_lemma/negpoint.py 1 1 0.5 1800 ;;
thmA)
  $PY apriori_regions/theorem_A/thmA_rig.py 1 1 0.5 0.796 7
  (wcopy apriori_regions/theorem_A && $PY fa_cert.py) ;;
lamstar)  # Theorem apriori: the 11 rows of lamstar.txt (runs_T0, T1, T2, B, B2), then the negative controls
  (wcopy apriori_regions/lamstar_sigma && rm -f runs_new.jsonl &&
   $PY batch.py 1.1 runs_new.jsonl 3.7498:4.0:14 4.0:4.45:14 4.45:5.05:13 5.05:5.8:12 5.8:6.75:11 6.75:8.05:10 \
       8.05:9.8:9 9.8:12.3:8 12.3:16.2:7 16.2:23.1:6 23.1:47.4:5 47.4:66:4 66:288:4 &&
   $PY negctl_row0.py) ;;
lambeta)  # Theorems aprioribeta and cont: recompute exactly the 120 + 25 rectangles of lamstar_beta.txt (-> certb/, certc/)
  (wcopy apriori_regions/lamstar_beta &&
   $PY -c "import json,glob; C=[json.load(open(f)) for f in sorted(glob.glob('$R/certificates/apriori_regions/apriori_cert/*.json'))]
json.dump([[d['s_lo'],d['s_hi'],d['beta_lo'],d['beta_hi'],d['Lstar']] for d in C if 'metodo' not in d], open('plan_from_certs.json','w'))
json.dump([[d['s_lo'],d['s_hi'],d['beta_lo'],d['beta_hi']] for d in C if 'metodo' in d], open('rects_from_certs.json','w'))" &&
   $PY run_beta.py plan_from_certs.json 0 1000 runsb_new.jsonl && $PY run_line.py rects_from_certs.json runs_line_new.jsonl &&
   $PY negctl_drv.py 30 31 0.35 0.4 4000; $PY negctl_drv.py 60 61 0 0.05 4000) ;;
route4)  # v1.1.0: joint constants of route 4 (66 certifications, ~45 CPU-min), T8 with sigma_min = 228, Theorem C on [228, 288]
  $PY g1_wall_lemma/route4/run_all.py
  $PY g0_constants/t8_sub.py 1 1 0.5 228 0 0.45
  G1_CRIT=crit4arb G1_MAXEVAL=3000 $PY g1_wall_lemma/driver.py 1 1 0.5 0.796 228 288 R4 1
  $PY g1_wall_lemma/negctl_crit4.py 1 1 0.5 225;  $PY g1_wall_lemma/negctl_crit4.py 1 1 0.5 224 0.3 0.45 ;;
g2)  # Appendix F: about 28.5 CPU-hours; on a cluster use g2_evans/slurm_array.sh instead (see README)
  (mkdir -p "$CAP_WORK/g2" && cp -R g2_evans/. "$CAP_WORK/g2/" && cd "$CAP_WORK/g2" && make evans &&
   for t in tasks_C.txt tasks_B.txt tasks_A.txt; do MODO=continuacao ./run_local.sh "$(nproc)" $t; done) ;;
*) echo "unknown step $s"; exit 2 ;;
esac; done
echo "done; outputs in $CAP_WORK"
