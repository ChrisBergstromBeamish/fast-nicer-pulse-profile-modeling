#!/usr/bin/env bash
# One step for the Desktop launcher: install the public code (if not there yet), then time it.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash "$HERE/install_reference_code.sh" | tee /tmp/physics_nicer_install_last.txt
if grep -q INSTALL_DONE /tmp/physics_nicer_install_last.txt; then
  bash "$HERE/run_public_code_timing.sh"
else
  echo "The install did not finish; the timing was not started."
fi
