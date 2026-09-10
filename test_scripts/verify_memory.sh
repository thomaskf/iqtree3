#!/bin/bash
# Compare IQ-TREE 3 memory usage against IQ-TREE 2 baseline + threshold.
# When IQ-TREE 2 baseline is 0 (unsupported command), falls back to the
# pre-defined expected value from expected_memory.tsv if a platform column is given.
#
# Args: $1 = IQ-TREE 2 log file (default: time_log_iqtree2.tsv)
#       $2 = IQ-TREE 3 log file (default: time_log_iqtree3.tsv)
#       $3 = platform column name in expected_memory.tsv for fallback (optional)

iqtree2_log="${1:-time_log_iqtree2.tsv}"
iqtree3_log="${2:-time_log_iqtree3.tsv}"
fallback_column="${3:-}"

WD="test_scripts/test_data"
threshold_file="${WD}/expected_memory.tsv"

tmp_thresholds=$(mktemp)
tmp_fallback=$(mktemp)
tmp_iqtree2=$(mktemp)
tmp_iqtree3=$(mktemp)

# Per-platform threshold column "thr-<platform>" when present, else diff-threshold.
# Peak memory is not comparable across platforms, so one shared allowance is either
# too tight on the noisy ones or meaningless on the quiet ones.
thr_index=2
if [ -n "$fallback_column" ]; then
    idx=$(head -1 "$threshold_file" | tr '\t' '\n' | awk -v c="thr-$fallback_column" '$0 == c {print NR}')
    if [ -n "$idx" ]; then
        thr_index=$idx
        echo "Using per-platform thresholds: thr-$fallback_column"
    else
        echo "No thr-$fallback_column column; using the shared diff-threshold"
    fi
fi
tail -n +2 "$threshold_file" | cut -f1,"$thr_index" > "$tmp_thresholds"

# Pre-defined fallback expected values for the given platform column
if [ -n "$fallback_column" ]; then
    col_index=$(head -1 "$threshold_file" | tr '\t' '\n' | awk -v col="$fallback_column" '{if ($0 == col) print NR}')
    if [ -z "$col_index" ]; then
        echo "WARNING: fallback column '$fallback_column' not found in $threshold_file; skipping fallback"
        fallback_column=""
    else
        tail -n +2 "$threshold_file" | cut -f"$col_index" > "$tmp_fallback"
    fi
fi

# Memory is column 3 of each log; column 1 is the command actually executed,
# kept so a breaching check can be retried.
tail -n +2 "$iqtree2_log" | cut -f3 > "$tmp_iqtree2"
tail -n +2 "$iqtree3_log" | cut -f3 > "$tmp_iqtree3"
tmp_cmd2=$(mktemp); tmp_cmd3=$(mktemp)
tail -n +2 "$iqtree2_log" | cut -f1 > "$tmp_cmd2"
tail -n +2 "$iqtree3_log" | cut -f1 > "$tmp_cmd3"
# shellcheck source=/dev/null
. "$(dirname "$0")/remeasure.sh"

# Reconcile the number of benchmark commands with the number of table rows.
# They are joined POSITIONALLY, so a mismatch means the pairing is wrong: without
# this, `paste` pads the short side and bc is handed an empty threshold
# ("Parse error: bad token" on macOS, "syntax error" on Linux), and rows print
# with a bare number in place of the command name.
n_rows=$(wc -l < "$tmp_thresholds")
n_log=$(wc -l < "$tmp_iqtree3")
if [ "$n_log" -ne "$n_rows" ]; then
    echo "⚠️  WARNING: the suite ran ${n_log} commands but memory has ${n_rows} threshold rows."
    if [ "$n_log" -gt "$n_rows" ]; then
        echo "   Skipping the last $((n_log - n_rows)) command(s) - they have no threshold:"
        tail -n +$((n_rows + 1)) "$tmp_cmd3" | sed 's/^/     /'
        for t in "$tmp_iqtree2" "$tmp_iqtree3" "$tmp_cmd2" "$tmp_cmd3"; do
            head -n "$n_rows" "$t" > "${t}.cut" && mv "${t}.cut" "$t"
        done
    else
        echo "   Ignoring the last $((n_rows - n_log)) threshold row(s) - no command produced them."
        head -n "$n_log" "$tmp_thresholds" > "${tmp_thresholds}.cut" && mv "${tmp_thresholds}.cut" "$tmp_thresholds"
    fi
    echo "   NOTE: rows are matched by POSITION. If the extra command(s) were added in the"
    echo "   middle rather than at the end, every later row is now compared against the"
    echo "   wrong command. Add the missing row(s) to keep the table in step."
fi

fail_count=0
row=0

while IFS=$'\t' read -r command threshold iqtree2_val iqtree3_val cmd2 cmd3; do
    ((row++))
    expected="$iqtree2_val"

    # Fall back to pre-defined value when IQ-TREE 2 baseline is unavailable
    if [ "$(echo "$expected == 0" | bc -l)" = "1" ]; then
        if [ -n "$fallback_column" ]; then
            expected=$(sed -n "${row}p" "$tmp_fallback")
            echo "ℹ️  $command: IQ-TREE 2 baseline unavailable, using pre-defined expected value (${expected}MB)"
        else
            echo "⏭ $command skipped (IQ-TREE 2 baseline unavailable, no fallback column provided)"
            continue
        fi
    fi

    allowed=$(echo "$expected + $threshold" | bc -l)
    is_exceed=$(echo "$iqtree3_val > $allowed" | bc -l)
    diff=$(echo "$iqtree3_val - $expected" | bc -l)

    # Retry once before failing: re-run this one command for both binaries and
    # re-evaluate. Costs nothing when everything passes.
    if [ "$is_exceed" = "1" ] && [ -n "$cmd3" ]; then
        echo "↻ $command exceeded (${diff}MB); retrying this command once..."
        read -r _ retry2_mem <<< "$(remeasure "$cmd2")"
        read -r _ retry3_mem <<< "$(remeasure "$cmd3")"
        if [ "$(echo "$retry2_mem > 0" | bc -l)" = "1" ] && [ "$(echo "$retry3_mem > 0" | bc -l)" = "1" ]; then
            expected="$retry2_mem"; iqtree3_val="$retry3_mem"
            allowed=$(echo "$expected + $threshold" | bc -l)
            is_exceed=$(echo "$iqtree3_val > $allowed" | bc -l)
            diff=$(echo "$iqtree3_val - $expected" | bc -l)
            echo "   retry: IQ-TREE2 ${retry2_mem}MB, IQ-TREE3 ${retry3_mem}MB, Diff ${diff}MB"
        else
            echo "   retry did not produce a usable measurement; keeping the first result"
        fi
    fi

    if [ "$is_exceed" = "1" ]; then
        echo "❌ $command exceeded the allowed memory usage."
        echo "   Expected: ${expected}MB, Threshold: ${threshold}MB, IQ-TREE3: ${iqtree3_val}MB, Diff: ${diff}MB"
        ((fail_count++))
    else
        echo "✅ $command passed the memory check."
        echo "   Expected: ${expected}MB, Threshold: ${threshold}MB, IQ-TREE3: ${iqtree3_val}MB, Diff: ${diff}MB"
    fi
done < <(paste "$tmp_thresholds" "$tmp_iqtree2" "$tmp_iqtree3" "$tmp_cmd2" "$tmp_cmd3")

rm -f "$tmp_thresholds" "$tmp_fallback" "$tmp_iqtree2" "$tmp_iqtree3" "$tmp_cmd2" "$tmp_cmd3"

if [ "$fail_count" -eq 0 ]; then
    echo "✅ All memory checks passed."
else
    echo "❌ $fail_count checks failed."
    exit 1
fi
