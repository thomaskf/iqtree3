#!/bin/bash
# =============================================================================
# benchmark_partitions.sh
#
# Run ModelFinder on each simulated partition dataset with different thread
# counts, record wall-clock time, and summarise results in a TSV file.
#
# Compares three IQ-TREE binaries:
#   IQTREE_CURRENT    – current (latest) modified version
#   IQTREE_FIRST_UPDATE – first update (intermediate modified version)
#   IQTREE_ORIG       – original unmodified IQ-TREE
#
# Usage:
#   bash benchmark_partitions.sh [IQTREE_CURRENT] [IQTREE_FIRST_UPDATE] [IQTREE_ORIG] [SIMDIR] [OUTDIR]
#
# Examples:
#   bash benchmark_partitions.sh ./iqtree3_current "" "" sim_partitions bench_partitions
#   bash benchmark_partitions.sh ./iqtree3_current ./iqtree3_first ./iqtree3_orig sim_partitions bench_partitions
# =============================================================================

IQTREE_CURRENT=${1:-iqtree3}
IQTREE_FIRST_UPDATE=${2:-}  # optional: first update (intermediate) binary
IQTREE_ORIG=${3:-}          # optional: original (unmodified) binary
SIMDIR=${4:-sim_partitions}
OUTDIR=${5:-bench_partitions}

THREAD_COUNTS=(1 2 4 8 12)

mkdir -p "$OUTDIR"

# -----------------------------------------------------------------------------
# Validate binaries
# -----------------------------------------------------------------------------
for bin in "$IQTREE_CURRENT" ${IQTREE_FIRST_UPDATE:+"$IQTREE_FIRST_UPDATE"} ${IQTREE_ORIG:+"$IQTREE_ORIG"}; do
    if [[ ! -x "$bin" ]] && ! command -v "$bin" &>/dev/null; then
        echo "ERROR: IQ-TREE binary not found or not executable: $bin"
        exit 1
    fi
done

# -----------------------------------------------------------------------------
# Load scenario metadata written by simulate_partitions.sh
# -----------------------------------------------------------------------------
META_TSV="$SIMDIR/scenarios.tsv"
if [[ ! -f "$META_TSV" ]]; then
    echo "ERROR: Scenario metadata not found: $META_TSV"
    echo "Run simulate_partitions.sh first."
    exit 1
fi

RESULTS="$OUTDIR/results.tsv"
echo -e "version\tscenario\tseq_type\tn_parts\ttotal_len\tmin_len\tmax_len\tnum_threads\twall_sec" > "$RESULTS"

echo "Current binary      : $IQTREE_CURRENT"
if [[ -n "$IQTREE_FIRST_UPDATE" ]]; then
    echo "First-update binary : $IQTREE_FIRST_UPDATE"
fi
if [[ -n "$IQTREE_ORIG" ]]; then
    echo "Original binary     : $IQTREE_ORIG"
fi
echo "Simulation dir  : $SIMDIR"
echo "Results         : $RESULTS"
echo ""

# -----------------------------------------------------------------------------
# Helper: time one ModelFinder run on a partitioned dataset.
# Args: label bin phy nex nt prefix scenario seq_type n_parts total_len min_len max_len
# -----------------------------------------------------------------------------
run_one() {
    local label="$1" bin="$2" phy="$3" nex="$4" nt="$5" prefix="$6"
    local scenario="$7" seq_type="$8" n_parts="$9" total_len="${10}" \
          min_len="${11}" max_len="${12}"

    # Pass sequence type explicitly so IQ-TREE interprets the alignment correctly
    local st_flag="-st $seq_type"

    printf "  [%-8s] -nt %2d ... " "$label" "$nt"

    local start end elapsed_ms elapsed_sec status
    start=$(python3 -c "import time; print(int(time.time() * 1000))")
    "$bin" -s "$phy" \
        -p "$nex" \
        -m TESTONLY \
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

    echo -e "$label\t$scenario\t$seq_type\t$n_parts\t$total_len\t$min_len\t$max_len\t$nt\t$elapsed_sec" >> "$RESULTS"
}

# -----------------------------------------------------------------------------
# Main benchmark loop — iterate over scenarios from metadata TSV
# -----------------------------------------------------------------------------
while IFS=$'\t' read -r scenario seq_type n_parts total_len min_len max_len phy nex; do
    # Skip header
    [[ "$scenario" == "scenario" ]] && continue
    # Skip failed simulations
    [[ "$phy" == "NA" ]] && { echo "  SKIP $scenario (simulation failed)"; continue; }
    # Skip if files are missing
    if [[ ! -f "$phy" || ! -f "$nex" ]]; then
        echo "  SKIP $scenario (files not found)"
        continue
    fi

    echo "--- $scenario  (type=$seq_type, $n_parts partitions, $total_len sites) ---"

    for nt in "${THREAD_COUNTS[@]}"; do
        prefix_cur="$OUTDIR/${scenario}_current_nt${nt}"
        run_one "current" "$IQTREE_CURRENT" "$phy" "$nex" "$nt" "$prefix_cur" \
            "$scenario" "$seq_type" "$n_parts" "$total_len" "$min_len" "$max_len"

        if [[ -n "$IQTREE_FIRST_UPDATE" ]]; then
            prefix_first="$OUTDIR/${scenario}_first_update_nt${nt}"
            run_one "first_update" "$IQTREE_FIRST_UPDATE" "$phy" "$nex" "$nt" "$prefix_first" \
                "$scenario" "$seq_type" "$n_parts" "$total_len" "$min_len" "$max_len"
        fi

        if [[ -n "$IQTREE_ORIG" ]]; then
            prefix_orig="$OUTDIR/${scenario}_original_nt${nt}"
            run_one "original" "$IQTREE_ORIG" "$phy" "$nex" "$nt" "$prefix_orig" \
                "$scenario" "$seq_type" "$n_parts" "$total_len" "$min_len" "$max_len"
        fi
    done
    echo ""
done < "$META_TSV"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo "============================================================"
echo "Benchmark complete."
echo "Results: $RESULTS"
echo ""
echo "Quick summary (wall time by scenario and thread count):"
echo "------------------------------------------------------------"
column -t -s $'\t' "$RESULTS"
