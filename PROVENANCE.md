# Provenance

All files were copied from the working tree of the project in which the proof was developed (October 2026). Apart
from the documentation, the scripts in `scripts/` and `g1_wall_lemma/negctl_crit3.py` (written for this
repository), every code, input and certificate file is a **byte-for-byte copy**, with one exception: 23 scripts
that contained absolute paths of the development machine. In those, a five-line header was inserted (it defines
the repository root, `certificates/` and the output directory `$CAP_WORK`, default `./work`), and each absolute
path was replaced by a path relative to the repository: inputs are read from `certificates/`, outputs are written
to `$CAP_WORK`. No other line was changed. One unit test was adapted: `g2_evans/tests/test_idempotente.py` used the
pilot certificates of G2, which are not distributed; it now uses `certificates/g2_evans/certs_final/` (one line).

| file | md5 of the original | md5 in this repository |
|---|---|---|
| `g0_sector_lemma/run_W.py` | `04f88bbec1d607fcb90ce14d8fdfc2d2` | `4e60cf8d1f28a2dabf6b7506d297d806` |
| `g0_sector_lemma/audit/negctl.py` | `2950fa7459904745d211785b1e6f5ad7` | `0bd16a8f17de1f7ab57e2783c7f0eade` |
| `g0_constants/cl.py` | `8d4d6e93955488d17b8c8fe8d7f85763` | `86e4de93b924d05a06fc8df534759fef` |
| `g0_constants/certline.py` | `170c2f05f9b5060b75073348339329ac` | `27630fb74eb4199ec94f2f4d230a9629` |
| `g0_constants/refine.py` | `d8c720c2595b4b8319b4c2ad955d8d30` | `1b6ae85ce68d9e0a2d5cd54bba1575a7` |
| `g0_constants/pbound.py` | `4a0cc171587ca830b6a7182890d4a9f5` | `cfa2b1f878a442d902ec04c95f1b2603` |
| `g0_constants/controls.py` | `150f6a68c9277257f8b7a8eedb0354a0` | `126fc6e534536f2c45894983bdb847b7` |
| `g0_constants/ptable.py` | `ea6796f22ded21df985600fd0d0594df` | `9dcde342cd111b1bd3d4db2659a123d4` |
| `g0_constants/tp_cert.py` | `b94bcf9dff11c8b0b7eb86ca2712c3a5` | `32b31d1521dcf97d4cbbfec97661736f` |
| `g0_constants/cert_T36.py` | `d460a5b59ff31615b610ef5d0c5539cc` | `e69144dad7cc67384487370608796536` |
| `g0_constants/cert_J2.py` | `23e0e0ee22d4e98a7c1ca2e1f32f18d9` | `743561616ba8aa70375a4216f3379e89` |
| `g0_constants/cert_T8.py` | `5d2dc0aa42eee632af07ccb8f1a168c6` | `fdc7e96127014a16cb281a21f26817c0` |
| `g0_constants/t8_patch.py` | `a22a3a921699b2dc4505e87145e4f527` | `a470b5e73c24a8657fe6ef721dae3d83` |
| `g0_constants/t8_sub.py` | `93459b423d8a9bbd7d872387df6e78b9` | `76ff251268b0a24095137a8c6b9d1436` |
| `g1_wall_lemma/arbcrit.py` | `650fbd4b5bf7630c240f47664a0d4742` | `ef8e6f8c17d2c12954b4c18a85ec0fb7` |
| `g1_wall_lemma/crit3arb.py` | `ac8cb586e9c496504361d54e72a1cfd0` | `725fb80da9176903837ddd1e90ea0c63` |
| `g1_wall_lemma/driver.py` | `3aebc605c53b2329e5a26d81445e33ef` | `89079423ff20653becb569a0e950a9f7` |
| `g1_wall_lemma/negpoint.py` | `e9eda49a6675c06279b3e4fa7212a6f0` | `0712dfdafdcd5696cd9db4a297e21704` |
| `g1_wall_lemma/heuristic/crit2.py` | `669d12b1797090fbc175323ec05742e3` | `bad13c12188e29972d132b55ab359501` |
| `g1_wall_lemma/heuristic/tp.py` | `3437f6746b354370630e39454cd42b18` | `66037c7a97338b43d7fbccd57bf4606d` |
| `g1_wall_lemma/pi_tails/cert_pi_arb.py` | `c3e4794cadec4eed1e9b04d121b90f52` | `3ee6d9ea71287cc4a5101cdcb0b98810` |
| `g1_wall_lemma/pi_tails/controles.py` | `91c4c1ce837492f93a7dc01024eeb392` | `e17a6c4118cfbfd2d889d709a15d0d13` |
| `g1_wall_lemma/pi_tails/audit_T36_caudas.py` | `7e22d576b92f9b47b74be75abeb644e4` | `6343756715637d556c825c497275d5dd` |

## Hashes quoted in the article (Remark 5.15)

| file | SHA-256 | here |
|---|---|---|
| `verify_final.txt` | `ba84b458863ccc7a3f2afdd5d62d4e435e6d0a5a492f36dc420c741c9d9f32a3` | `certificates/g2_evans/verify_final.txt` |
| `g2_certificados_final.tar.gz` | `0b3fde2f7cb27b4773e1ea7466ee8ea1a8f3096a71232b02548b56b80429e4aa` | not distributed; member-by-member manifest in `certificates/g2_evans/ARCHIVE_MANIFEST.txt` |
| `lamstar.txt` | `77907aa1deadf37a…` | `certificates/apriori_regions/lamstar.txt`, `g2_evans/lamstar.txt` |
| `lamstar_beta.txt` | `b8a3e384f93bb87f…` | `certificates/apriori_regions/lamstar_beta.txt`, `g2_evans/lamstar_beta.txt` |
| `sweep.py` | `920145f9…` | `g2_evans/sweep.py` |
| `evans.c` | `176d5693…` | `g2_evans/src/evans.c` |

All are checked by `scripts/check_certificates.py`. The archive `g2_certificados_final.tar.gz` (189 files: the G2
certificates, `verify_final.txt`, the inputs, `sweep.py` and `src/` of the final run, and development copies of five G1
scripts) is **not** included: those five scripts contain absolute paths of the development machine (four of them)
or differ in an output path (t8_sub.py). `ARCHIVE_MANIFEST.txt` lists the SHA-256 of every member; the other 184
members are in this repository byte for byte, which `check_certificates.py` verifies. The repository archive
(`micropolar-couette-cap-v1.0.0.zip`) is a different file with its own checksum.

## Code versions recorded in the certificates

* The 145 rectangle certificates of `lamstar_beta.txt` record the md5 of the code that produced them
  (`ap_arb.py` 5db5a481…, `ap_airy.py` e846c156…, `ap_airy3.py` 9d30867d…, `driver2.py` 77ab5d71…,
  `driver3.py` c2ac7137…, `driver5.py` 78c4bf3f…). These are exactly the files in `apriori_regions/lamstar_beta/`.
* `lamstar.txt` (runs `runs_T0/T1/T2/B/B2.jsonl`) was produced with an earlier `ap_arb.py` (md5 653ed993…) run
  through `batch.py` (which disables the probing heuristic). That version is kept in `apriori_regions/lamstar_sigma/`
  under its original module name. Recomputing the first row with it reproduces `runs_T0.jsonl` exactly
  (18 160 evaluations, 9 363 boxes); running `driver2.py` directly (probing on) also certifies the row, with other counts.
* Theorem C: the runs `out_C3R2`, `out_C3A2` (code v3), `out_C4R2`, `out_C4R2b` (v4) and the bands `out_S_*` (v5) were
  made with successive versions of `crit3arb.py`/`driver.py`; only the last version is kept here. Recomputed here: `out_C4R2b` is reproduced exactly; `out_C3A2` is reproduced with the same
  342 boxes, 347 evaluations and 0 failures, but a larger minimal margin (0.433 instead of 0.028), because the
  parameter heuristic now also offers the relative-cone variants of Lemma F8; `out_S_280` is reproduced with 0 failures
  and 65 boxes instead of 67, because the T8 table with σ_min = 215 was added after that run.
* `T8_cert_1_1_0.5_s215.json` was produced by two passes of `t8_sub.py`; the tolerance of the second pass (cells
  [1.7, 2.2] and [2.2, 3]) was not logged.
* `certificates/g1_wall_lemma/tail_3e5.log` and `certificates/apriori_regions/thmA_cA_alpha7.log` were regenerated for
  this repository: the original runs printed to the terminal only.

## Software

FLINT 3.0.1 (Ubuntu 24.04 package `libflint-dev 3.0.1-3.1build1`), GMP 6.3.0, MPFR 4.2.1, gcc 13.3.0 for `evans`;
Python 3.13.16, python-flint 0.9.0 (its wheel bundles FLINT 3.6.0), numpy 2.5.3, scipy 1.18.1, mpmath 1.3.0.

## Not included (not on the logical path of the Main Theorem)

Theorem W and the cone rows η = 0, 0.75, 0.885–0.8895; Romanov's constant a(0); the uniformity statements in g and N;
the other triples (N,J,G); Lemma F7 (route 4, Lemma F5, was added in v1.1.0, see below); everything labelled *numerical* (sections on numerics and
on instability for G > J); the floating-point predecessors of the verifiers; the pilot runs of G2; positive controls
and cross-validations that need the floating-point collocation solver. These are in the project archive, not here.

## v1.1.0 (route 4)

All files added in v1.1.0 were written for this version (none is a copy of a development file), except
`g1_wall_lemma/heuristic/Mtab_1_1_0.5.json` (floating-point table of the project, md5 of the original kept unchanged).
The route-4 certificates were produced on 2026-10-29 with the code in `g1_wall_lemma/route4/` and `g1_wall_lemma/crit4arb.py`
(md5 recorded in `M_cert_1_1_0.5.json`), python-flint 0.9.0, numpy 2.5.3.  After the runs, only the module docstring
of `mjoint_cert.py` was corrected (grid (Z/128)^2 instead of (Z/64)^2, usage line); the md5 in `M_cert_1_1_0.5.json` is
that of the corrected file (the runs used md5 13ffa15a…), and the sample M_u(0) was recomputed with it (identical).  The relative-cone table
`T8_cert_1_1_0.5_s228.json` was produced by the unchanged `g0_constants/t8_sub.py 1 1 0.5 228 0 0.45`.  No file of
v1.0.0 was modified except the documentation, `scripts/` and `CHANGELOG.md`.
