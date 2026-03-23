#!/bin/bash
# =============================================================================
# benchmark_modelfinder.sh
#
# Run ModelFinder on each simulated alignment with different thread counts,
# record wall-clock time, and summarise the results in a TSV file.
#
# Optionally accepts a second IQ-TREE binary (unmodified/current version) to
# run a side-by-side comparison against the modified version.
#
# Usage:
#   bash benchmark_modelfinder.sh [IQTREE_NEW] [SIMDIR] [OUTDIR] [IQTREE_ORIG]
#
# Examples:
#   bash benchmark_modelfinder.sh ./iqtree3 ./sim_data ./bench_results
#   bash benchmark_modelfinder.sh ./iqtree3_new ./sim_data ./bench_results ./iqtree3_orig
# =============================================================================

IQTREE_NEW=${1:-iqtree3}
SIMDIR=${2:-sim_alignments}
OUTDIR=${3:-bench_results}
IQTREE_ORIG=${4:-}          # optional: original (unmodified) binary

# Thread counts to test
THREAD_COUNTS=(1 2 4 8 16)

mkdir -p "$OUTDIR"

# Validate binaries
for bin in "$IQTREE_NEW" ${IQTREE_ORIG:+"$IQTREE_ORIG"}; do
    if [[ ! -x "$bin" ]] && ! command -v "$bin" &>/dev/null; then
        echo "ERROR: IQ-TREE binary not found or not executable: $bin"
        exit 1
    fi
done

RESULTS="$OUTDIR/results.tsv"

# Write TSV header — includes 'version' column to distinguish modified vs original
echo -e "version\talignment\tdata_type\tnsites\tnum_threads\teffective_threads_expected\twall_sec" > "$RESULTS"

echo "Modified binary : $IQTREE_NEW"
if [[ -n "$IQTREE_ORIG" ]]; then
    echo "Original binary : $IQTREE_ORIG"
fi
echo "Simulation dir  : $SIMDIR"
echo "Results         : $RESULTS"
echo ""

# Helper: infer data type and number of states from filename
# Prints: <seqtype>  <nstate>
parse_alignment() {
    local base
    base=$(basename "$1" .phy)
    if [[ "$base" == dna_* ]];   then echo "DNA 4";    return; fi
    if [[ "$base" == aa_* ]];    then echo "AA 20";    return; fi
    if [[ "$base" == codon_* ]]; then echo "CODON 61"; return; fi
    echo "AUTO 0"
}

# Helper: read number of sites from a PHYLIP file (first line, second field)
get_nsites() {
    awk 'NR==1{print $2}' "$1" 2>/dev/null
}

# Helper: time a single ModelFinder run and append one row to $RESULTS
# Args: version_label  iqtree_bin  aln  seqtype  nsites  nstate  nt  expected  prefix
run_one() {
    local label="$1" bin="$2" aln="$3" seqtype="$4"
    local nsites="$5" nstate="$6" nt="$7" expected="$8" prefix="$9"

    local st_flag=""
    [[ "$seqtype" != "AUTO" ]] && st_flag="-st $seqtype"

    echo -n "  [${label}] -nt $nt  (expected effective: $expected) ... "

    local start end elapsed_ms elapsed_sec status
    start=$(python3 -c "import time; print(int(time.time() * 1000))")
    "$bin" -s "$aln" \
        -m TEST \
        -nt "$nt" \
        $st_flag \
        --prefix "$prefix" \
        --redo -quiet 2>/dev/null
    status=$?
    end=$(python3 -c "import time; print(int(time.time() * 1000))")

    elapsed_ms=$(( end - start ))
    elapsed_sec=$(python3 -c "print(f'{$elapsed_ms / 1000:.3f}')")

    if [[ $status -eq 0 ]]; then
        echo "${elapsed_sec}s"
    else
        echo "FAILED"
        elapsed_sec="NA"
    fi

    local base
    base=$(basename "$aln" .phy)
    echo -e "$label\t$base\t$seqtype\t$nsites\t$nt\t$expected\t$elapsed_sec" >> "$RESULTS"
}

# =============================================================================
# Main benchmark loop
# =============================================================================
for aln in "$SIMDIR"/*.phy; do
    [[ -f "$aln" ]] || continue

    read -r seqtype nstate <<< "$(parse_alignment "$aln")"
    nsites=$(get_nsites "$aln")
    base=$(basename "$aln" .phy)

    echo "--- $base  (type=$seqtype, nsites=$nsites) ---"

    for nt in "${THREAD_COUNTS[@]}"; do
        # Expected effective threads after cap: max(1, nptn*nstate/4000), capped at nt
        local_expected=$nt
        if [[ "$nstate" -gt 0 ]]; then
            local_expected=$(( nsites * nstate / 4000 ))
            (( local_expected < 1 )) && local_expected=1
            (( local_expected > nt )) && local_expected=$nt
        fi

        run_one "modified" "$IQTREE_NEW" "$aln" "$seqtype" \
            "$nsites" "$nstate" "$nt" "$local_expected" \
            "$OUTDIR/${base}_modified_nt${nt}"

        if [[ -n "$IQTREE_ORIG" ]]; then
            # Original always uses all nt threads (no cap), so expected == nt
            run_one "original" "$IQTREE_ORIG" "$aln" "$seqtype" \
                "$nsites" "$nstate" "$nt" "$nt" \
                "$OUTDIR/${base}_original_nt${nt}"
        fi
    done
    echo ""
done

echo "============================================================"
echo "Benchmark complete."
echo "Results: $RESULTS"
echo ""
echo "Quick summary (wall time by alignment and thread count):"
echo "------------------------------------------------------------"
column -t -s $'\t' "$RESULTS"
