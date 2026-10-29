# Negative controls

A negative control runs a verifier on a target that is slightly **false** (or just outside the certified range) and
checks that it is **not** certified. Logs written in Portuguese use *NAO certifica / falha* = not certified / fails.

| component | target that must fail | where | result |
|---|---|---|---|
| Sector lemma (App. A) | θ₂ = 0.53803 (true max 0.5380375), m = 0.6338 (true min 0.633763), and for η = 0.75 θ₂ = 0.8052, θ₁ = −0.0934 | `g0_sector_lemma/audit/negctl.py` → `g0/negctl_sector_lemma.log` | 4/4 not certified (irreducible box < 1e-6 at the extremum) |
| Boundary constants (App. B) | ã(0.05) > 0.35346272 | `g0_constants/refine.py 0.05 ut 0.35346272` → `constants/refine_at_0.05_neg.log` | not certified |
| Cone at η = 0.05 / 0.1 (App. B) | θ₂ = 0.54935 / 0.56135 | `g0_constants/controls.py` (needs scipy) | not certified (as recorded in the project reports; rerun with `scripts/full_recompute.sh consts`) |
| Π constants (App. D) | Π = 0.1308120942 (below a rigorous point value); tail lemmas with T = 2, 5, 10 | `g1_wall_lemma/pi_tails/controles.py` → `../g1_wall_lemma/pi_tails/controles.log` | not certifiable / lemmas inapplicable |
| Theorem C, refined criterion (App. D) | σ = 260 and σ = 275, β ∈ [0, 1e-9], all 32 parameter choices | `g1_wall_lemma/negctl_crit3.py` → `g1/negctl_theoremC.log` | not certified |
| Theorem C with route 4 (App. D, v1.1.0) | σ = 225, β ∈ [0, 1e-9]; σ = 224, β ∈ [0.2, 0.3] and [0.3, 0.45], all 32 parameter choices (positive control: σ = 228 certified, margin 0.109) | `g1_wall_lemma/negctl_crit4.py` → `route4/negctl_route4.log` | not certified |
| Route-4 joint constants (App. D, v1.1.0) | M_u(0) ≤ 0.2313 | lower enclosure 0.23138 in `certificates/g1_wall_lemma/route4/M_cert_1_1_0.5.json` (rigorous grid value) | refuted |
| Λ*(σ) (App. E) | Λ* − 1 in every row (point scan); full verifier with Λ* = 13 on [3.7498, 4] | `apriori_regions/lamstar_sigma/negctl.py`, `negctl_row0.py` → `lamstar_sigma/` | fails in all 11 rows; 29 729 failing boxes |
| Λ*(σ,β) and empty region (App. E) | β 0.1–0.2 below each certified limit (point criterion on the line); full verifier on [30,31]×[0.35,0.4] and [60,61]×[0,0.05] | `apriori_regions/lamstar_beta/negctl_beta.py`, `negctl_cont2.py`, `negctl_drv.py` → `lamstar_beta/` | fails (ρ ≥ 1); 242 and 14 failing boxes |
| Evans function (App. F) | unit tests of the anchoring logic: path through a Theorem-A box, corner contact only, component of line boxes only | `g2_evans/tests/test_ancoragem.py` (`make test`) | rejected |

**Known discrepancy with the article.** Appendix D says that at β ∈ [0, 10⁻⁹] no parameter choice certifies σ = 280
or 260. With the code in this repository (and with the code of the original run), σ = 260 and σ = 275 fail, but
at σ = 280 the choice λ = 1, θ_ρ = 0.8, variant (iir) gives a positive rigorous **point** margin 0.0116 (also at
σ = 278: 0.0066). This does not affect the proof (σ₀′ = 288 is certified on boxes, and [215, 288] is only used through
the β-bands), but the sentence in the article should read "260 or 275" instead of "280 or 260". (Corrected in the article after v1.0.0.)
