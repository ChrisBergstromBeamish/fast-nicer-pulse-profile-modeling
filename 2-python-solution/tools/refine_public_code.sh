#!/usr/bin/env bash
# The public code X-PSI with its settings refined ONE AT A TIME from its shipped ones: its energies first, then its rotation
# steps (leaves), then its mesh (cells), then its rays; then two finer combinations. For both shipped checks, one
# evaluation at each setting. Its value at the check's star and its bins are saved (results/public_code_refine.jsonl in
# 2-python-solution, and an .npz file of bins beside it for each setting); the times it writes are from one evaluation
# each and are not used. compare_with_public_code.py then sets ours beside them.
# Needs install_reference_code.sh to have ended with INSTALL_DONE. Written and run on 23 Nov 2026, after the first run of
# the public code on this machine.
set -u
BASE=/root/physics_nicer
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$HERE/../results/public_code_refine.jsonl"
export PATH="$BASE/miniforge/bin:$PATH"
export OMP_NUM_THREADS=1
mkdir -p "$HERE/../results"
mkdir -p "$BASE/run" && cd "$BASE/run" || { echo "REFINE_FAILED no run folder"; exit 20; }
{
echo "[$(date -u +%H:%M:%S)] start: settings are sqrt_num_cells, leaves, rays, energies"
for S in "32 64 512 128 1 quick" "32 64 512 192 1 quick" "32 64 512 256 1 quick" "32 64 512 384 1 quick" "32 64 512 512 1 quick" \
         "32 128 512 64 1 quick" "32 256 512 64 1 quick" \
         "64 64 512 64 1 quick" "128 64 512 64 1 quick" \
         "32 64 1024 64 1 quick" "32 64 2048 64 1 quick" \
         "64 128 512 384 1 quick" "128 256 512 384 1 quick" \
         "32 64 512 192 1 second" "32 64 512 256 1 second" "32 64 512 384 1 second" "32 64 512 512 1 second" \
         "32 128 512 128 1 second" "32 256 512 128 1 second" \
         "64 64 512 128 1 second" "128 64 512 128 1 second" \
         "32 64 1024 128 1 second" "32 64 2048 128 1 second" \
         "64 128 512 384 1 second" "128 256 512 384 1 second"; do
  echo "[$(date -u +%H:%M:%S)] settings: $S"
  conda run --no-capture-output -n xpsi python "$HERE/time_public_code.py" "$BASE/xpsi" "$OUT" $S 2>&1 | grep -E "check ;|^at p:|TIMING_DONE|Error|Traceback|rror:" || echo "REFINE_FAILED at settings $S"
done
echo "[$(date -u +%H:%M:%S)] end"
echo ALL_REFINE_DONE
} 2>&1 | tee -a "$BASE/refine.log"
