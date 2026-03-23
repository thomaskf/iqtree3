#!/bin/bash
# =============================================================================
# benchmark_modelfinder.sh
#
# Run ModelFinder on each simulated alignment with different thread counts,
# record wall-clock time, and summarise the results in a TSV file.
#
# Usage:
#   bash benchmark_modelfinder.sh [IQTREE_BIN] [SIMDIR] [OUTDIR]
#
# Examples:
#   bash benchmark_modelfinder.sh
#   bash benchmark_modelfinder.sh /usr/local/bin/iqtree3 ./sim_alignments ./bench_results
# =============================================================================

IQTREE=${1:-iqtree3}
SIMDIR=${2:-sim_alignments}
OUTDIR=${3:-bench_results}
RESULTS="$OUTDIR/results.tsv"

# Thread counts to test
THREAD_COUNTS=(1 2 4 8 16)

mkdir -p "$OUTDIR"

if ! command -v "$IQTREE" &>/dev/null; then
    echo "ERROR: IQ-TREE binary not found: $IQTREE"
    exit 1
fi

# Write TSV header
echo -e "alignment\tdata_type\tnsites\tnum_threads\teffective_threads_expected\twall_sec" > "$RESULTS"

echo "Results will be written to: $RESULTS"
echo ""

# Helper: infer data type and expected effective thread cap from filename
# Returns: <seqtype_flag>  <nstate>
parse_alignment() {
    local f="$1"
    local base
    base=$(basename "$f" .phy)

    if [[ "$base" == dna_* ]];   then echo "DNA 4";   return; fi
    if [[ "$base" == aa_* ]];    then echo "AA 20";   return; fi
    if [[ "$base" == codon_* ]]; then echo "CODON 60"; return; fi
    echo "AUTO 0"
}

# Helper: read number of sites from a PHYLIP file (first line, second field)
get_nsites() {
    awk 'NR==1{print $2}' "$1" 2>/dev/null
}

# Main benchmark loop
for aln in "$SIMDIR"/*.phy; do
    [[ -f "$aln" ]] || continue

    read -r seqtype nstate <<< "$(parse_alignment "$aln")"
    nsites=$(get_nsites "$aln")
    base=$(basename "$aln" .phy)

    # Estimate number of patterns (upper bound = nsites; in practice fewer)
    # Use nsites as a proxy since we don't know nptn without loading the file
    nptn_approx=$nsites

    echo "--- $base  (type=$seqtype, nsites=$nsites) ---"

    for nt in "${THREAD_COUNTS[@]}"; do
        # Expected effective threads after cap: max(1, nptn*nstate/4000), capped at nt
        if [[ "$nstate" -gt 0 ]]; then
            expected=$(( nptn_approx * nstate / 4000 ))
            (( expected < 1 )) && expected=1
            (( expected > nt )) && expected=$nt
        else
            expected=$nt
        fi

        prefix="$OUTDIR/${base}_nt${nt}"

        st_flag=""
        if [[ "$seqtype" != "AUTO" ]]; then
            st_flag="-st $seqtype"
        fi

        echo -n "  -nt $nt  (expected effective: $expected) ... "

        # Time the ModelFinder run
        start=$(date +%s%3N)
        "$IQTREE" -s "$aln" \
            -m TEST \
            -nt "$nt" \
            $st_flag \
            --prefix "$prefix" \
            --redo -quiet 2>/dev/null
        status=$?
        end=$(date +%s%3N)

        elapsed_ms=$(( end - start ))
        elapsed_sec=$(echo "scale=3; $elapsed_ms / 1000" | bc)

        if [[ $status -eq 0 ]]; then
            echo "${elapsed_sec}s"
        else
            echo "FAILED"
            elapsed_sec="NA"
        fi

        echo -e "$base\t$seqtype\t$nsites\t$nt\t$expected\t$elapsed_sec" >> "$RESULTS"
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
