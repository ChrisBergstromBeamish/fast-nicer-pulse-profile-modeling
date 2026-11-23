#!/usr/bin/env bash
# Times the compiled program on the public quick example at several settings (one process, one thread).
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
echo "machine:$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) ; compiler: $(gcc --version | head -n 1)"
for s in "12 16 24" "10 16 24" "10 16 16" "8 12 16" "24 24 32"; do
  echo "settings (rings per panel, points, energies): $s"
  ./nicer_fast tables.bin $s 400 1.4 12 1.0 0.5000000000000001 0.0 1.2217304763960306 0.75 6.7
done
echo "a different star at every working-out (10 rings per panel, 16 points, 24 energies; then 8, 12, 16):"
./nicer_fast tables.bin 10 16 24 -400 1.4 12 1.0 0.5000000000000001 0.0 1.2217304763960306 0.75 6.7
./nicer_fast tables.bin 8 12 16 -400 1.4 12 1.0 0.5000000000000001 0.0 1.2217304763960306 0.75 6.7
