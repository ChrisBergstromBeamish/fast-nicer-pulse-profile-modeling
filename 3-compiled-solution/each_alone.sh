#!/usr/bin/env bash
# Each star of a stars file worked out alone: one run of the program per star (check.py compares this with the batch mode).
# usage, in this folder on the Linux side: bash each_alone.sh <program> <rings per panel> <points> <energies> <stars file> <out stem>
# <stem>.txt gets one line per star: every bit of ln L (a hexadecimal floating-point number) and the code 0, or "refused" and the
# refusal's code; <stem>.bin gets the dumps (bin counts, then model rates) of the stars that are not refused, one after another.
prog=$1; npp=$2; n=$3; nE=$4; stars=$5; out=$6
dump="$out.one.bin"
: > "$out.txt"; : > "$out.bin"
while read -r line; do
  [ -z "$line" ] && continue
  res=$(./"$prog" tables.bin "$npp" "$n" "$nE" 0 $line "$dump")
  case "$res" in
    "refused: code "*) echo "refused ${res#refused: code }" >> "$out.txt" ;;
    *) echo "$(echo "$res" | sed -n 's/^every bit of ln L: //p') 0" >> "$out.txt"; cat "$dump" >> "$out.bin" ;;
  esac
  rm -f "$dump"
done < "$stars"
