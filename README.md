# micropolar-couette-cap

Codes and certificates for the computer-assisted proof of the **Main Theorem** of

> *Linear stability of plane Couette flow of a micropolar fluid with the microrotation–vorticity boundary condition*
> P. Braz e Silva, preprint, 2026.

**Main Theorem (computer-assisted).** Let $(N,J,G)=(1,1,\tfrac12)$. Consider two-dimensional perturbations of micropolar
plane Couette flow with $w=\tfrac12\,\mathrm{curl}\,u$ on the walls. Then for all $\alpha>0$ and $R>0$ every eigenvalue
satisfies $\alpha R\,C_i\le-\tfrac1{20}$. In particular, the flow is spectrally stable for every Reynolds number.

The proof combines analytic arguments, given in the article, with finite verifications in ball arithmetic
(Arb/FLINT). This repository contains **exactly the computational links of the logical path** of the proof
(article, Remark 5.16 and Table "The logical path of the Main Theorem"), their certificates and their negative
controls. It contains nothing else: no prototypes, no other triples $(N,J,G)$, no numerical evidence.

## Map: ingredient of the logical path → code → certificate

| ingredient (article) | code | certificate(s) | time to recompute (2 cores) |
|---|---|---|---|
| Sector lemma at η = 0.05 (and 0.1): $\mathcal L\ne0$, cone $\cos\theta\ge0.545117$, $\lvert\mathcal G\rvert\ge0.6083$ (App. A) | `g0_sector_lemma/run_W.py`, `g0lib.py` | `certificates/g0_sector_lemma/out_C1_005.json`, `out_C1_010.json` | 25 s + 27 s |
| $a(0.05)$, $\tilde a(0.05)$, $m_F$, $P(S')$ (App. B) | `g0_constants/certline.py`, `refine.py`, `pbound.py`, `cl.py` | `certificates/g0_constants/line_0.05.json`, `log_P_0.05.txt` | ≈ 2 min |
| Tables T1–T2 ($\Omega_A,\Omega_\psi$), T3–T6 ($K_\Pi,K_{\Pi1},Q,Q^*$), T7 (fine $P(S')$), T8 (relative cone, Lemma F8), $J_G,J_H,J_{sH}$ (App. D) | `g0_constants/tp_cert.py`, `cert_T36.py`, `ptable.py`, `cert_T8.py`, `t8_patch.py`, `t8_sub.py`, `cert_J2.py` | `certificates/g0_constants/tp_cert_0.05.json`, `T36_cert_0.05.json`, `ptable_0.05.json`, `T8_cert_1_1_0.5.json`, `T8_cert_1_1_0.5_s215.json`, `J_cert_1_1_0.5.json` | ≈ 30 min (T8: 15 min) |
| $\Pi,\Pi',\Pi_\psi,\Pi'_\psi$ with tails in Arb (App. D, Lemma C2) | `g1_wall_lemma/pi_tails/cert_pi_arb.py` | `certificates/g1_wall_lemma/cert_pi_arb.json` (+ the values used: `cert_pi_0.05.json`, `cert_pipsi.json`) | seconds |
| $\mathfrak K>0$ on $\sigma\in[288,3\cdot10^5]$, $0\le\beta\le\beta_A$; tail lemma $\sigma\ge3\cdot10^5$ (Theorem C, App. D) | `g1_wall_lemma/driver.py` with `crit3arb.py` / `arbcrit.py`; `tail.py` | `certificates/g1_wall_lemma/out_C4R2b.json`, `out_C4R2.json`, `out_C3R2.json`, `out_C3A2.json`, `out_A2.json` | ≈ 13 min |
| Theorem C in β-bands, σ ∈ [215, 288] (only used to discard boxes of K′) | same, `G1_BMIN` | `certificates/g1_wall_lemma/strips/out_S_*.json`, `g1strip.txt` | ≈ 17 min |
| $c_A=0.796$ for α ≥ 7; pairs $F_A\ge\sigma_k$ on $[\alpha_k,7]$; $F_A(0^+)\in[3.749882242,3.749882243]$ (Remark cA) | `apriori_regions/theorem_A/thmA_rig.py`, `fa_cert.py` | `certificates/apriori_regions/fa_cert.json` (c_A: output of `thmA_rig.py`, see below) | seconds |
| A priori region $\Lambda^*(\sigma)$ (Theorem apriori, App. E) | `apriori_regions/lamstar_sigma/batch.py`, `driver2.py`, `ap_arb.py` | `certificates/apriori_regions/lamstar.txt`, `lamstar_runs/runs_*.jsonl` | ≈ 1 h |
| A priori region $\Lambda^*(\sigma,\beta)$ and empty region (Theorem aprioribeta); exact Airy criterion on the line + continuation (Theorem cont) (App. E) | `apriori_regions/lamstar_beta/run_beta.py`, `driver3.py`, `run_line.py`, `driver5.py`, `ap_arb.py`, `ap_airy*.py` | `certificates/apriori_regions/lamstar_beta.txt`, `apriori_cert/*.json` (145) | ≈ 4 h |
| Computation G2: validated Evans function on $\mathcal K'\cap\{\sigma<85,\ \beta<\beta'_\varnothing\}$ (App. F; Remark 5.15) | `g2_evans/src/evans.c` (C, FLINT), `g2_evans/sweep.py` | `certificates/g2_evans/certs_final/*.jsonl` (609 records, 523 valid), `verify_final.txt`, `ARCHIVE_MANIFEST.txt` | ≈ 28.5 CPU-hours |

Negative controls (targets that must **not** be certified) are listed in [`certificates/negative_controls/README.md`](certificates/negative_controls/README.md).
Everything that is proved analytically (model, Theorem A, the analytic part of Theorem C, the Evans function and the
continuation principle) is in the article only.

## Requirements

* **Level 1 (quick verify):** Python 3, standard library only (tested with Python 3.13). No compiler, no FLINT. Runs in a few seconds.
* **Level 2/3 (recomputation):** FLINT **3.0.1** (C library, used by `g2_evans/evans`), GMP, MPFR, a C compiler;
  Python 3 with **python-flint 0.9.0**, numpy, scipy, mpmath (the original runs used Python 3.13.16, numpy 2.5.3,
  scipy 1.18.1, mpmath 1.3.0, gcc 13.3.0 on Ubuntu 24.04). Note that the python-flint 0.9.0 wheel bundles its own
  FLINT (3.6.0); the Python certificates were produced with that wheel. See `environment.yml` and `Dockerfile`.

## Reproduction

### Level 1 — quick verify (seconds, standard library only)

```sh
sh scripts/quick_verify.sh
```

It (1) checks the SHA-256 of every code, input and certificate file against `SHA256SUMS` (documentation and
metadata files — `*.md`, `LICENSE*`, `CITATION.cff`, `.zenodo.json`, `Dockerfile`, `environment.yml` — are not
hashed, so that they can be completed without regenerating it); (2) rereads every certificate of the logical path and
checks the claims made in the article (`scripts/check_certificates.py`); (3) runs the final G2 verification
`sweep.py verify` on the recorded Evans certificates and compares its output line by line with
`certificates/g2_evans/verify_final.txt`. The last line must be

```
QUICK VERIFY PASSED: K' COBERTO, 0 autovalores instaveis (K' covered, 0 unstable eigenvalues)
```

(`sweep.py` prints its report in Portuguese; *K′ coberto* = K′ covered, *autovalores instáveis* = unstable eigenvalues.)
Level 1 checks the certificates **as recorded**. It does not redo interval computations.

### Level 2 — recompute a sample of every component (≈ 5 minutes on 2 cores)

```sh
python3 scripts/recompute_samples.py      # outputs go to ./work (or $CAP_WORK); certificates/ is never modified
```

Clean build of `evans`, three G2 boxes (two line, one contour), the sector lemma at η = 0.05 (full run), $P(S')$, T1–T6, $J$,
the bound for $\tilde a(0.05)$, three full runs of the Theorem C criterion, the tail lemma, $\Pi$ constants with tails,
$c_A$, $F_A$, one row of $\Lambda^*(\sigma)$, one rectangle of $\Lambda^*(\sigma,\beta)$, one continuation rectangle, and
negative controls. Every recorded field except wall-clock times is compared with the certificate. `evans` is built
with the original flags (`-O2 -march=native`, x86-64 with FMA); on other hardware the floating-point bookkeeping may
subdivide differently, which changes counts but not the validity of a certificate (that is what `sweep.py verify` checks).

### Level 3 — full recomputation

`scripts/full_recompute.sh` lists and runs every command of the logical path, component by component
(`sh scripts/full_recompute.sh g0 consts pi g1 thmA lamstar lambeta g2`). Approximate CPU times: G0 1 min, constants
30 min, Π seconds, G1 30 min, Theorem A seconds, Λ*(σ) 1 h, Λ*(σ,β) 4 h, **G2 28.5 h**.

**G2 on a cluster (SLURM).** In `g2_evans/`: `make` (set `FLINT_PREFIX` if FLINT is not in `/usr`), then for each task
list in the order C → B → A

```sh
sbatch --array=0-49 --export=ALL,BOXES=tasks_C.txt,NTASKS=50,MODO=continuacao slurm_array.sh
```

The sweep is resumable (rerun the same array after a time limit). Then
`python3 sweep.py verify --boxes boxes.txt --certs 'certs/*.jsonl'`. Without SLURM: `MODO=continuacao ./run_local.sh <cores> tasks_C.txt`.

### Docker

```sh
docker build --target verify -t cap-verify . && docker run --rm cap-verify     # level 1
docker build --target full   -t cap-full   . && docker run --rm cap-full       # level 2 (FLINT 3.0.1 + python-flint 0.9.0)
```

## Layout

```
g0_sector_lemma/    Appendix A: sector lemma (integral representation, tails, affine box envelopes)
g0_constants/       Appendices B, D: a, a~, m_F, P(S'), tables T1-T8, J (all on the line z_eta, eta = 0.05)
g1_wall_lemma/      Appendix D: criterion of Theorem C in Arb (crit3arb, arbcrit), driver, tail lemma
  heuristic/          floating-point code that only CHOOSES free parameters of the criterion (not part of the rigour)
  pi_tails/           Pi constants with tails in Arb
apriori_regions/    Remark cA and Appendix E: Theorem A constants, Lambda*(sigma), Lambda*(sigma,beta), exact-Airy continuation
g2_evans/           Appendix F: evans.c (FLINT), sweep.py, base boxes and task lists, and the inputs exactly as used in the final run
certificates/       all final certificates, by component (+ negative_controls/)
tests/              run_tests.sh (level 1; --flint: unit tests of evans/sweep; --samples: level 2)
scripts/            quick_verify.sh, check_certificates.py, verify_sha256sums.py, recompute_samples.py, full_recompute.sh
PROVENANCE.md       origin of every file, code hashes, the path-only edits, and known differences
```

## Caveats (article, Remark 5.17)

1. The verification codes have been tested (positive and negative controls, comparisons with independent
   floating-point solvers, unit tests of the anchoring logic) but have **not** had an independent line-by-line review.
2. The proof relies on the correctness of FLINT/Arb and of the compiler. The a priori verifiers also use cellwise
   IEEE double arrays with outward rounding implemented by hand; the geometric bookkeeping of `sweep.py` uses floating
   point with relative tolerance 1e-12 on boxes of a single canonical subdivision tree.
3. The analytic proofs have been checked only within the project and have not been refereed.
4. `sweep.py verify` (level 1) checks the certificates as recorded; it does not re-evaluate the Evans function.
   Re-running G2 takes about 28 CPU-hours.

See also `PROVENANCE.md` for the few places where the current code reproduces a certificate without matching it
byte for byte (the certificate is reproduced, i.e. recomputed without failure; counts differ because the code gained
options after the run).

## How to cite

Please cite the article and this archive (see `CITATION.cff`; the DOI will be assigned by Zenodo):

> P. Braz e Silva. *micropolar-couette-cap: codes and certificates for the computer-assisted proof of linear stability
> of micropolar plane Couette flow*, version 1.0.0, 2026. Zenodo. doi:10.5281/zenodo.TODO

## License

Code: MIT (`LICENSE`). Certificates and data (`certificates/`, `g2_evans/*.txt`, `g2_evans/apriori_cert/`,
`g2_evans/g1cert/`, `g2_evans/fa_cert.json`, `g1_wall_lemma/heuristic/*.json|*.npz`): CC BY 4.0 (`LICENSE-data`).
