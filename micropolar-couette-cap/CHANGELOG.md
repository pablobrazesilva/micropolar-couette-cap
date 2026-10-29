# Changelog

## v1.1.0 — 2026-10-29

Route 4 certified (optional strengthening; the Main Theorem and its logical path are unchanged).

* `g1_wall_lemma/route4/`: certification of the joint constants of Lemma F5 for (1,1,1/2) — `mjoint_cert.py` (sup over
  the plane (t_k, t_f) of |f_beta|, by a strip with second-order cells, integration-by-parts tails in t_f, a uniform
  bound in t_f and a quasi-static far field in t_k), `linetab.py`, `ftab.py` (Arb tables), `run_all.py`.
* `g1_wall_lemma/crit4arb.py`: criterion of Theorem C with routes 3 and 4 (crit3arb itself is unchanged), and an extra
  relative-cone table with sigma_min = 228; `g1_wall_lemma/negctl_crit4.py`; `heuristic/Mtab_1_1_0.5.json` (floating-point
  joint constants, used only to choose the free parameters).
* `certificates/g1_wall_lemma/route4/`: `M_cert_1_1_0.5.json` (M_u, M_v on 32 nodes beta in [0,1], M^s_u, M^s_v, with
  two-sided enclosures), the 66 run records, `T8_cert_1_1_0.5_s228.json`, and `out_R4.json`: K > 0 on
  sigma in [228, 288], 0 <= beta <= beta_A (309 boxes, 0 failures).  New certified threshold sigma_0' = 228
  (alpha_0' = 16.93) instead of 288.
* `certificates/negative_controls/route4/`.
* `scripts/check_certificates.py`, `scripts/recompute_samples.py` (component `route4`), `scripts/full_recompute.sh`
  (step `route4`) extended; `PROVENANCE.md` updated.

## v1.0.0 — 2026-10-28

First public release, accompanying the preprint.

* Codes and certificates of every computational link of the logical path of the Main Theorem for (N,J,G) = (1,1,1/2):
  sector lemma (App. A), boundary constants and tables T1–T8, J (App. B, D), Π constants with tails in Arb,
  criterion of Theorem C on σ ∈ [288, 3·10⁵] plus tail lemma and the β-bands on [215, 288] (App. D), Theorem A
  constants, a priori regions Λ*(σ), Λ*(σ,β), the empty region and the exact-Airy continuation (App. E), and the
  validated Evans-function sweep of K′ (App. F; package v6, 523 certificates).
* Negative controls for every component.
* `scripts/quick_verify.sh` (standard library only), `scripts/recompute_samples.py`, `scripts/full_recompute.sh`,
  Dockerfile with targets `verify` and `full`.
