#!/usr/bin/env bash
# Times the batch mode of nicer_fast (the fast form) on many stars at once, on 1, 2, 4, 6 and 8 threads. In this folder, on the
# Linux side, after make:   bash time_batch.sh
# 1. 4,096 stars drawn from the public example's prior (draw_stars.py, seed 1), at the settings (10 rings per panel, 16 points,
#    24 energies) and (8, 12, 16): five repeats per thread count; stars per second, the fastest and the median of the five, the
#    speed-up over one thread and that speed-up divided by the number of threads.
# 2. The batch mode on one thread against the program's own one-at-a-time timing (a negative number of repeats): the same 2,800
#    stars (the 7 blocks of 400 that mode works out, written to a file by the awk below), time per star; then, in five rounds,
#    the two with the same statistic (the fastest of the same 7 blocks of 400) one right after the other.
# 3. What limits the speed-up: one thread held on each processor in turn, and four threads held on processors 0 to 3 and then 4
#    to 7 (taskset), at (10, 16, 24), on the first 1,024 of the 4,096 stars.
# 4. The same question without threads: 4 and then 8 separate programs at the same moment, each on one thread held on its own
#    processor (0 to 3, then 0 to 7), each working out the same 1,024 stars; the sum of their stars per second.
# usage: bash time_batch.sh [parts]      (the parts to run, for example "4"; all four if none is named)
# The power state is read from the Windows side before and after each part, and after each timed line (2 = mains, 1 = battery;
# "unknown" off Windows).
# The time of a batch is the program's own (the stars only: not the start of the program or the loading of the tables).
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
PYTHON3=${PYTHON3:-python3}
REPEATS=5
power() { local s; s=$(powershell.exe -NoProfile -Command "(Get-CimInstance Win32_Battery).BatteryStatus" 2>/dev/null | tr -d '\r'); echo "power state ${s:-unknown} ($(date '+%Y-%m-%d %H:%M:%S %Z'))"; }
state() { local s; s=$(powershell.exe -NoProfile -Command "(Get-CimInstance Win32_Battery).BatteryStatus" 2>/dev/null | tr -d '\r'); echo "power ${s:-unknown} $(date '+%H:%M:%S')"; }
# five repeats of one batch: prints "fastest median" in stars per second; the arguments are the command before the stars file
repeat_batch() {
  local stars=$1; shift
  for r in $(seq $REPEATS); do "$@" "$stars" 2>&1 >/dev/null | sed -n 's/.* \([0-9.]*\) stars per second$/\1/p'; done |
    sort -g | awk '{v[NR] = $1} END {printf "%.1f %.1f\n", v[NR], (NR % 2 ? v[(NR + 1) / 2] : 0.5 * (v[NR / 2] + v[NR / 2 + 1]))}'
}
PARTS=${*:-1 2 3 4}
part() { case " $PARTS " in *" $1 "*) return 0 ;; esac; return 1; }
echo "machine:$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) ; processors: $(nproc) ; compiler: $(gcc --version | head -n 1)"
[ -x nicer_fast ] && [ -f tables.bin ] || { echo "run make first"; exit 2; }
"$PYTHON3" -I -B draw_stars.py 1 4096 stars_timing.txt
head -n 1024 stars_timing.txt > stars_timing_1024.txt
if part 1; then
echo "== 1. 4096 stars, the fast form, $REPEATS repeats each: stars per second, fastest and median ; speed-up (by the fastest) ; speed-up / threads"
power
for s in "10 16 24" "8 12 16"; do
  ./nicer_fast tables.bin batch $s 1 stars_timing.txt > batch_timing_out.txt 2>/dev/null
  echo "setting (rings per panel, points, energies) = ($s): refused $(awk '$2 != 0' batch_timing_out.txt | wc -l) of $(wc -l < batch_timing_out.txt), early exit (ln L = -inf) $(awk '$1 == "-inf" && $2 == 0' batch_timing_out.txt | wc -l)"
  one=""
  for th in 1 2 4 6 8; do
    read -r fast med < <(repeat_batch stars_timing.txt ./nicer_fast tables.bin batch $s $th)
    [ -z "$one" ] && one=$fast
    echo "  threads $th: $fast fastest, $med median ; speed-up $(awk -v a=$fast -v b=$one 'BEGIN {printf "%.2f", a / b}') ; per thread $(awk -v a=$fast -v b=$one -v t=$th 'BEGIN {printf "%.2f", a / b / t}') ; $(state)"
  done
done
power
fi
if part 2; then
echo "== 2. the batch mode on one thread against the one-at-a-time timing (the same 2800 stars), time per star"
power
for s in "10 16 24" "8 12 16"; do
  awk 'BEGIN {split("1.4 12 1.0 0.5000000000000001 0.0 1.2217304763960306 0.75 6.7", p, " ");
              for (k = 0; k < 7; k++) for (i = 0; i < 400; i++) { line = "";
                for (m = 0; m < 8; m++) { q = (m == 4) ? 1.0e-3 * sin(0.9 * i + k) : p[m + 1] * (1.0 + 1.0e-3 * sin(1.7 * i + 2.3 * m + 0.4 * k));
                  line = line sprintf("%.17g ", q) } print line } }' > stars_one_at_a_time.txt
  read -r fast med < <(repeat_batch stars_one_at_a_time.txt ./nicer_fast tables.bin batch $s 1)
  echo "  ($s) batch on one thread, 2800 stars, $REPEATS repeats: $(awk -v f=$fast 'BEGIN {printf "%.4f", 1000 / f}') ms per star fastest, $(awk -v f=$med 'BEGIN {printf "%.4f", 1000 / f}') ms median ; $(state)"
  for r in $(seq $REPEATS); do
    echo "  ($s) one at a time: $(./nicer_fast tables.bin $s -400 1.4 12 1.0 0.5000000000000001 0.0 1.2217304763960306 0.75 6.7 | grep 'time per' | cut -d';' -f1); $(state)"
  done
  # The same statistic for both, at nearly the same moment: in each round, the one-at-a-time timing (the fastest of its 7
  # blocks of 400), then the batch on one thread on each of the same 7 blocks of 400, the fastest of them.
  awk '{print > ("stars_one_at_a_time_block" int((NR - 1) / 400) ".txt")}' stars_one_at_a_time.txt
  for r in $(seq $REPEATS); do
    one=$(./nicer_fast tables.bin $s -400 1.4 12 1.0 0.5000000000000001 0.0 1.2217304763960306 0.75 6.7 | sed -n 's/.*each time: \([0-9.]*\) ms.*/\1/p')
    bat=$(for k in 0 1 2 3 4 5 6; do ./nicer_fast tables.bin batch $s 1 stars_one_at_a_time_block$k.txt 2>&1 >/dev/null | sed -n 's/.* \([0-9.]*\) stars per second$/\1/p'; done |
          sort -g | tail -n 1 | awk '{printf "%.4f", 1000 / $1}')
    echo "  ($s) round $r, the fastest of the same 7 blocks of 400: one at a time $one ms, batch on one thread $bat ms, batch minus one at a time $(awk -v a=$bat -v b=$one 'BEGIN {printf "%+.1f%%", 100 * (a - b) / b}') ; $(state)"
  done
  rm -f stars_one_at_a_time_block*.txt
done
power
fi
if part 3; then
echo "== 3. what limits the speed-up, at (10, 16, 24), the first 1024 of the stars, $REPEATS repeats: stars per second, fastest and median"
power
for cpu in 0 1 2 3 4 5 6 7; do
  read -r fast med < <(repeat_batch stars_timing_1024.txt taskset -c $cpu ./nicer_fast tables.bin batch 10 16 24 1)
  echo "  one thread held on processor $cpu: $fast fastest, $med median ; $(state)"
done
for cpus in 0-3 4-7; do
  read -r fast med < <(repeat_batch stars_timing_1024.txt taskset -c $cpus ./nicer_fast tables.bin batch 10 16 24 4)
  echo "  four threads held on processors $cpus: $fast fastest, $med median ; $(state)"
done
power
fi
if part 4; then
echo "== 4. separate programs at the same moment, one thread each held on its own processor, at (10, 16, 24), the first 1024 stars each, $REPEATS repeats: the sum of their stars per second, fastest and median"
power
for n in 4 8; do
  for r in $(seq $REPEATS); do
    for cpu in $(seq 0 $((n - 1))); do
      taskset -c $cpu ./nicer_fast tables.bin batch 10 16 24 1 stars_timing_1024.txt 2> "batch_timing_p$cpu.txt" > /dev/null &
    done
    wait
    cat batch_timing_p*.txt | sed -n 's/.* \([0-9.]*\) stars per second$/\1/p' | awk '{s += $1} END {printf "%.1f\n", s}'
    rm -f batch_timing_p*.txt
  done | sort -g | awk -v n=$n '{v[NR] = $1} END {printf "  %d programs on processors 0 to %d: %.1f fastest, %.1f median ; ", n, n - 1, v[NR], v[(NR + 1) / 2]}'
  echo "$(state)"
done
power
fi
rm -f batch_timing_out.txt stars_one_at_a_time.txt stars_timing_1024.txt
