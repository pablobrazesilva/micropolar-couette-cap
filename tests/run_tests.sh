#!/bin/sh
# Test suite of the repository.
#   level 1 (standard library only):  sh tests/run_tests.sh
#   + FLINT unit tests of evans:      sh tests/run_tests.sh --flint
#   + sample recomputation (level 2): sh tests/run_tests.sh --flint --samples
set -e
cd "$(dirname "$0")/.."
PY="${PYTHON:-python3}"
sh scripts/quick_verify.sh
case " $* " in *" --flint "*)
  # unit tests of the Evans code and of sweep.py: reference values, idempotent resume, anchoring logic (V2),
  # including the negative cases (path through a Theorem-A box, corner contact only, component without anchor)
  make -C g2_evans test ;;
esac
case " $* " in *" --samples "*) "$PY" scripts/recompute_samples.py ;; esac
echo "TESTS PASSED"
