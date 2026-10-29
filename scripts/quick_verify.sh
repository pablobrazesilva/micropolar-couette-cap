#!/bin/sh
# Level 1 (quick) verification: Python >= 3.8 standard library only, no compilation, a few minutes.
#  1. SHA-256 of every code, input and certificate file (SHA256SUMS);
#  2. the recorded certificates of every component of the logical path (scripts/check_certificates.py);
#  3. the final G2 verification `sweep.py verify` on the recorded Evans certificates, compared line by line with
#     certificates/g2_evans/verify_final.txt (the output quoted in the article).
# Nothing is written to disk.
set -e
PY0="${PYTHON:-python3}"; case "$PY0" in */*) PY0="$(cd "$(dirname "$PY0")" && pwd)/$(basename "$PY0")";; esac   # absolute path, resolved before cd
cd "$(dirname "$0")/.."
PY="$PY0"
echo "== 1/3 file integrity"
"$PY" -I scripts/verify_sha256sums.py
echo "== 2/3 component certificates"
"$PY" -I scripts/check_certificates.py
echo "== 3/3 G2: sweep.py verify on certificates/g2_evans/certs_final (re-checks every certificate, (V1), (V2), coverage of K')"
cd g2_evans
OUT="$("$PY" -I sweep.py verify --boxes boxes.txt --certs '../certificates/g2_evans/certs_final/*.jsonl')"
echo "$OUT"
if [ "$OUT" = "$(cat ../certificates/g2_evans/verify_final.txt)" ]; then
  echo "== output identical to certificates/g2_evans/verify_final.txt"
else
  echo "== OUTPUT DIFFERS from certificates/g2_evans/verify_final.txt"; exit 1
fi
echo "$OUT" | grep -q "^RESULTADO: K' COBERTO, 0 autovalores instaveis$"
echo "QUICK VERIFY PASSED: K' COBERTO, 0 autovalores instaveis (K' covered, 0 unstable eigenvalues)"
