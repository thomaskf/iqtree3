#!/bin/bash
# =============================================================================
# benchmark_modelfinder.sh
#
# Run ModelFinder on each simulated alignment with different thread counts,
# comparing the original IQ-TREE against the updated version at each value
# of the --mf-thread-factor denominator j.
#
# Usage:
#   bash benchmark_modelfinder.sh [IQTREE_UPDATED] [IQTREE_ORIG] [SIMDIR] [OUTDIR]
#
# Examples:
#   bash benchmark_modelfinder.sh ./iqtree3_updated ./iqtree3_orig sim_alignments bench_mf
# =============================================================================

IQTREE_UPDATED=${1:-iqtree3}
IQTREE_ORIG=${2:-}          # optional: original (unmodified) binary
SIMDIR=${3:-sim_alignments}
OUTDIR=${4:-bench_mf}

# Thread counts to test
THREAD_COUNTS=(1 2 4 8 12)

# Values of j (denominator in max(1, nptn*nstate/j)) to test for updated binary
J_VALUES=(1000 2000 3000 4000 5000 6000 7000)

mkdir -p "$OUTDIR"

# -----------------------------------------------------------------------------
# Validate binaries
# -----------------------------------------------------------------------------
for bin in "$IQTREE_UPDATED" ${IQTREE_ORIG:+"$IQTREE_ORIG"}; do
    if [[ ! -x "$bin" ]] && ! command -v "$bin" &>/dev/null; then
        echo "ERROR: IQ-TREE binary not found or not executable: $bin"
        exit 1
    fi
done

RESULTS="$OUTDIR/results.tsv"
echo -e "version\talignment\tdata_type\tnsites\tj_factor\tnum_threads\twall_sec" > "$RESULTS"

echo "Updated binary : $IQTREE_UPDATED"
if [[ -n "$IQTREE_ORIG" ]]; then
    echo "Original binary: $IQTREE_ORIG"
fi
echo "Simulation dir : $SIMDIR"
echo "Results        : $RESULTS"
echo ""

# -----------------------------------------------------------------------------
# Helper: infer data type from filename prefix
# -----------------------------------------------------------------------------
parse_seqtype() {
    local base
    base=$(basename "$1" .phy)
    if [[ "$base" == dna_* ]];   then echo "DNA";   return; fi
    if [[ "$base" == aa_* ]];    then echo "AA";    return; fi
    if [[ "$base" == codon_* ]]; then echo "CODON"; return; fi
    echo "AUTO"
}

# -----------------------------------------------------------------------------
# Helper: read number of sites from a PHYLIP file (first line, second field)
# -----------------------------------------------------------------------------
get_nsites() {
    awk 'NR==1{print $2}' "$1" 2>/dev/null
}

# -----------------------------------------------------------------------------
# Helper: time a single ModelFinder run and append one row to $RESULTS
# Args: version_label  iqtree_bin  aln  seqtype  nsites  nt  j_factor  prefix
# (j_factor = "NA" for original which has no --mf-thread-factor support)
# -----------------------------------------------------------------------------
run_one() {
    local label="$1" bin="$2" aln="$3" seqtype="$4"
    local nsites="$5" nt="$6" j_factor="$7" prefix="$8"

    local extra_flags=""
    [[ "$j_factor" != "NA" ]] && extra_flags="--mf-thread-factor $j_factor"

    printf "  [%-14s] -nt %2d  j=%-5s ... " "$label" "$nt" "$j_factor"

    local start end elapsed_ms elapsed_sec status
    start=$(python3 -c "import time; print(int(time.time() * 1000))")
    "$bin" -s "$aln" \
        -m TESTONLY \
        -nt "$nt" \
        -st "$seqtype" \
        $extra_flags \
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
    echo -e "$label\t$base\t$seqtype\t$nsites\t$j_factor\t$nt\t$elapsed_sec" >> "$RESULTS"
}

# =============================================================================
# Main benchmark loop
# =============================================================================
for aln in "$SIMDIR"/*.phy; do
    [[ -f "$aln" ]] || continue

    seqtype=$(parse_seqtype "$aln")
    nsites=$(get_nsites "$aln")
    base=$(basename "$aln" .phy)

    echo "--- $base  (type=$seqtype, nsites=$nsites) ---"

    for nt in "${THREAD_COUNTS[@]}"; do
        # Original binary: no j factor (label as NA)
        if [[ -n "$IQTREE_ORIG" ]]; then
            run_one "original" "$IQTREE_ORIG" "$aln" "$seqtype" \
                "$nsites" "$nt" "NA" \
                "$OUTDIR/${base}_original_nt${nt}"
        fi

        # Updated binary: test each j value
        for j in "${J_VALUES[@]}"; do
            run_one "updated_j${j}" "$IQTREE_UPDATED" "$aln" "$seqtype" \
                "$nsites" "$nt" "$j" \
                "$OUTDIR/${base}_updated_j${j}_nt${nt}"
        done
    done
    echo ""
done

echo "============================================================"
echo "Benchmark complete."
echo "Results: $RESULTS"
echo ""
echo "Quick summary:"
echo "------------------------------------------------------------"
column -t -s $'\t' "$RESULTS"
