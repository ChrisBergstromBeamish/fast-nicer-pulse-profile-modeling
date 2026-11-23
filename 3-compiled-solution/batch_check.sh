#!/usr/bin/env bash
# The batch mode next to each star alone, for check.py (which compares what this writes). In this folder, on the Linux side:
#     bash batch_check.sh <program> <rings per panel> <points> <energies> <stem> [<rings per panel> <points> <energies> <stem> ...]
# For each group: <stem>.txt is worked out alone (each_alone.sh) and in batches on 1, 2, 3, 4 and 8 threads, <stem>_shuffled.txt
# (the same stars in another order) in batches on 1, 2, 3, 4 and 8 threads, and <stem>.txt once more on 8 threads. Everything
# goes to batch_check_out/<program>_<rings>_<points>_<energies>_<stem>_*.
prog=$1; shift
mkdir -p batch_check_out
while [ $# -ge 4 ]; do
  npp=$1; n=$2; nE=$3; stem=$4; shift 4
  o="batch_check_out/${prog}_${npp}_${n}_${nE}_${stem}"
  bash each_alone.sh "$prog" "$npp" "$n" "$nE" "$stem.txt" "${o}_alone"
  for th in 1 2 3 4 8; do
    ./"$prog" tables.bin batch "$npp" "$n" "$nE" "$th" "$stem.txt" "${o}_b$th.bin" > "${o}_b$th.txt" 2>/dev/null
    ./"$prog" tables.bin batch "$npp" "$n" "$nE" "$th" "${stem}_shuffled.txt" "${o}_s$th.bin" > "${o}_s$th.txt" 2>/dev/null
  done
  ./"$prog" tables.bin batch "$npp" "$n" "$nE" 8 "$stem.txt" "${o}_again.bin" > "${o}_again.txt" 2>/dev/null
done
