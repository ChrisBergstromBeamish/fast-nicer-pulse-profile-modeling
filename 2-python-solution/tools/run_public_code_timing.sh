#!/usr/bin/env bash
# Times the public code X-PSI on its own quick example, at its shipped settings and at finer ones (to see where its
# value goes as its settings are refined). Needs install_reference_code.sh to have ended with INSTALL_DONE.
# NOT RUN AND NOT TESTED when written (22 Nov 2026).
set -u
BASE=/root/physics_nicer
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$HERE/../results/public_code_times.jsonl"
export PATH="$BASE/miniforge/bin:$PATH"
export OMP_NUM_THREADS=1
mkdir -p "$BASE/run" && cd "$BASE/run" || { echo "TIMING_FAILED no run folder"; exit 20; }
{
echo "[$(date -u +%H:%M:%S)] machine: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) ; $(nproc) threads"
echo "[$(date -u +%H:%M:%S)] the public code's own tests of its two shipped likelihood values and of its likelihood step"
( cd "$BASE/run" && conda run --no-capture-output -n xpsi python -m pytest -q "$BASE/xpsi/xpsi/tests/test_likelihood.py" "$BASE/xpsi/xpsi/tests/test_background_marginalisation.py" ) || echo "OWN_TESTS_FAILED"
for S in "32 64 512 64 12 quick" "64 128 512 128 6 quick" "128 256 512 128 3 quick" "192 384 512 192 1 quick"          "32 64 512 128 4 second" "64 128 512 128 2 second" "128 256 512 128 1 second" "192 384 512 192 1 second"; do
  echo "[$(date -u +%H:%M:%S)] settings: $S"
  conda run --no-capture-output -n xpsi python "$HERE/time_public_code.py" "$BASE/xpsi" "$OUT" $S || echo "TIMING_FAILED at settings $S"
done
echo "[$(date -u +%H:%M:%S)] end"
echo ALL_TIMINGS_DONE
} 2>&1 | tee -a "$BASE/timing.log"
