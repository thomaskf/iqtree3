#!/bin/bash
# =============================================================================
# simulate_partitions.sh
#
# Simulate multi-partition alignments using IQ-TREE AliSim for benchmarking
# the new thread-scheduling strategy in PartitionFinder.
#
# Each scenario varies the number of partitions and their sizes (nptn * nstate)
# to exercise different regions of the scheduler's behaviour:
#
#   Scenario A – Many uniform light partitions   (16 x 500 bp DNA)
#   Scenario B – Mixed heavy + light partitions  (1 x 8000 bp + 15 x 200 bp DNA)
#   Scenario C – Few heavy partitions            (4 x 4000 bp DNA)
#   Scenario D – Many heavy partitions           (8 x 2000 bp DNA)
#   Scenario E – Very many very light partitions (32 x 100 bp DNA)
#
# Usage:
#   bash simulate_partitions.sh [IQTREE_BIN] [OUTDIR]
#
# Examples:
#   bash simulate_partitions.sh
#   bash simulate_partitions.sh ./iqtree3 sim_partitions
# =============================================================================

IQTREE=${1:-iqtree3}
OUTDIR=${2:-sim_partitions}
NTAXA=20
SEED=12345
DNA_MODEL="GTR+G"

mkdir -p "$OUTDIR"

# Verify IQ-TREE binary
if ! command -v "$IQTREE" &>/dev/null && [[ ! -x "$IQTREE" ]]; then
    echo "ERROR: IQ-TREE binary not found: $IQTREE"
    exit 1
fi

echo "IQ-TREE binary : $IQTREE"
echo "Output dir     : $OUTDIR"
echo "Taxa: $NTAXA  Seed: $SEED"
echo ""

# -----------------------------------------------------------------------------
# Helper: write a NEXUS partition file.
# Usage: write_nexus <nexus_path> <model> <len1> <len2> ...
# -----------------------------------------------------------------------------
write_nexus() {
    local nex_path="$1"
    local model="$2"
    shift 2
    local lengths=("$@")
    local n=${#lengths[@]}
    local pos=1

    {
        echo "#nexus"
        echo "begin sets;"
        for ((i = 0; i < n; i++)); do
            local end=$(( pos + lengths[i] - 1 ))
            printf "  charset part%d = %d-%d;\n" $((i + 1)) "$pos" "$end"
            pos=$(( end + 1 ))
        done
        # Build charpartition line
        local cp=""
        for ((i = 0; i < n; i++)); do
            [[ $i -gt 0 ]] && cp+=", "
            cp+="${model}:part$((i + 1))"
        done
        echo "  charpartition myparts = $cp;"
        echo "end;"
    } > "$nex_path"
}

# -----------------------------------------------------------------------------
# Helper: run AliSim for one scenario and print result.
# Usage: run_alisim <scenario_name> <len1> <len2> ...
# Returns: 0 on success, 1 on failure
# Sets global PHY_OUT and NEX_OUT to the output file paths.
# -----------------------------------------------------------------------------
run_alisim() {
    local sname="$1"
    shift
    local lengths=("$@")
    local n_parts=${#lengths[@]}
    local total_len=0
    for l in "${lengths[@]}"; do (( total_len += l )); done
    local min_len=${lengths[0]}
    local max_len=${lengths[0]}
    for l in "${lengths[@]}"; do
        (( l < min_len )) && min_len=$l
        (( l > max_len )) && max_len=$l
    done

    local prefix="$OUTDIR/$sname"
    NEX_OUT="${prefix}_partition.nex"
    PHY_OUT="${prefix}.phy"

    write_nexus "$NEX_OUT" "$DNA_MODEL" "${lengths[@]}"

    "$IQTREE" --alisim "$prefix" \
        -p "$NEX_OUT" \
        -t "RANDOM{yh/$NTAXA}" \
        --seed "$SEED" \
        --redo -quiet 2>/dev/null

    if [[ $? -eq 0 && -f "$PHY_OUT" ]]; then
        printf "  OK   %-25s  %2d partitions  total %6d bp\n" \
            "$sname" "$n_parts" "$total_len"
        # Append metadata row: scenario n_parts total_len min_len max_len phy nex
        printf "%s\t%d\t%d\t%d\t%d\t%s\t%s\n" \
            "$sname" "$n_parts" "$total_len" "$min_len" "$max_len" \
            "$PHY_OUT" "$NEX_OUT" >> "$META_TSV"
        return 0
    else
        printf "  FAIL %-25s\n" "$sname"
        printf "%s\t%d\t%d\t%d\t%d\tNA\tNA\n" \
            "$sname" "$n_parts" "$total_len" "$min_len" "$max_len" >> "$META_TSV"
        return 1
    fi
}

# -----------------------------------------------------------------------------
# Write metadata TSV header
# -----------------------------------------------------------------------------
META_TSV="$OUTDIR/scenarios.tsv"
echo -e "scenario\tn_parts\ttotal_len\tmin_len\tmax_len\tphy\tnex" > "$META_TSV"

# -----------------------------------------------------------------------------
# Define and run scenarios
# -----------------------------------------------------------------------------
echo "=== Simulating partition alignments ==="

# Scenario A: 16 uniform light partitions (16 × 500 bp)
A_LENS=()
for ((i = 0; i < 16; i++)); do A_LENS+=(500); done
run_alisim "A_uniform_light" "${A_LENS[@]}"

# Scenario B: Mixed (1 × 8000 bp + 15 × 200 bp)
B_LENS=(8000)
for ((i = 0; i < 15; i++)); do B_LENS+=(200); done
run_alisim "B_mixed" "${B_LENS[@]}"

# Scenario C: Few heavy partitions (4 × 4000 bp)
C_LENS=()
for ((i = 0; i < 4; i++)); do C_LENS+=(4000); done
run_alisim "C_few_heavy" "${C_LENS[@]}"

# Scenario D: Many heavy partitions (8 × 2000 bp)
D_LENS=()
for ((i = 0; i < 8; i++)); do D_LENS+=(2000); done
run_alisim "D_many_heavy" "${D_LENS[@]}"

# Scenario E: Very many very light partitions (32 × 100 bp)
E_LENS=()
for ((i = 0; i < 32; i++)); do E_LENS+=(100); done
run_alisim "E_many_very_light" "${E_LENS[@]}"

echo ""
echo "Metadata saved to: $META_TSV"
echo ""
echo "Scenario summary:"
column -t -s $'\t' "$META_TSV"
