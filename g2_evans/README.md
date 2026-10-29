# g2_evans — validated Evans function and the sweep of K′ (Appendix F)

This directory is the package of the final G2 run (version v6), with its inputs exactly as used.

* `src/evans.c`, `src/tm.inc`, `src/lohner.inc` — Evans function of the 6th-order system through the compound matrix
  Λ³C⁶ (dimension 20), Taylor integrator with a majorant tail bound, Taylor models in (Λ, β², S) with Lohner/QR
  remainder control, argument principle on B(Λ*) (contour mode) or certification of the line Re Λ = −δ′_inf
  (line mode). Built with `make` (FLINT ≥ 3.0, `FLINT_PREFIX=/path` if not in /usr); `make test` runs the unit tests.
* `sweep.py` — base boxes of K′ (`boxes`), adaptive runs (`run --modo contorno|continuacao`) and the final check
  (`verify`): validity of every rectangle certificate of `lamstar_beta.txt`, of every G1 band of `g1strip.txt`, of
  every Evans certificate, conditions (V1)–(V2) by union–find over face adjacency, and the coverage of every base box
  through the canonical subdivision tree. Standard library only.
* Inputs: `boxes.txt` (263 base boxes), `tasks_A/B/C.txt`, `lamstar.txt`, `lamstar_beta.txt` + `apriori_cert/`,
  `g1strip.txt` + `g1cert/`, `fa_cert.json` (identical copies of the files in `../certificates/`).
* `slurm_array.sh`, `run_local.sh` — cluster / local drivers (environment variable `MODO=continuacao` for line mode).

The certificates of the final run are in `../certificates/g2_evans/certs_final/`; to check them:

    python3 sweep.py verify --boxes boxes.txt --certs '../certificates/g2_evans/certs_final/*.jsonl'

which must reproduce `../certificates/g2_evans/verify_final.txt` line by line. Messages are in Portuguese
(*K' COBERTO, 0 autovalores instaveis* = K′ covered, 0 unstable eigenvalues).
