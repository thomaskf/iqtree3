#!/bin/bash
# =============================================================================
# simulate_partitions.sh
#
# Simulate multi-partition alignments using IQ-TREE AliSim for benchmarking
# the partition thread-scheduling strategy.
#
# Three data types (DNA, AA, CODON) × six scenarios per type = 18 datasets.
#
# Partition sizes (same for all data types; divisible by 3 for CODON):
#   very_light = 102   sites (codons for CODON)
#   very_heavy = 10002 sites (codons for CODON)
#
# Scenarios (a-f):
#   a – 2  very light partitions
#   b – 10 very light partitions
#   c – 2  very heavy partitions
#   d – 10 very heavy partitions
#   e – 2  mixed (1 very heavy + 1 very light)
#   f – 10 mixed (5 very heavy + 5 very light)
#
# Usage:
#   bash simulate_partitions.sh [IQTREE_BIN] [OUTDIR]
# =============================================================================

IQTREE=${1:-iqtree3}
OUTDIR=${2:-sim_partitions}
NTAXA=20
SEED=12345

mkdir -p "$OUTDIR"

if ! command -v "$IQTREE" &>/dev/null && [[ ! -x "$IQTREE" ]]; then
    echo "ERROR: IQ-TREE binary not found: $IQTREE"
    exit 1
fi

echo "IQ-TREE binary : $IQTREE"
echo "Output dir     : $OUTDIR"
echo "Taxa: $NTAXA  Seed: $SEED"
echo ""

# -----------------------------------------------------------------------------
# Data types: parallel arrays (index-aligned)
# -----------------------------------------------------------------------------
DTYPE_NAMES=("DNA"   "AA"   "CODON")
DTYPE_MODELS=("GTR+G" "LG+G" "GY+G")

# Partition sizes — same for all data types; divisible by 3 so CODON works
VERY_LIGHT=102
VERY_HEAVY=10002

# -----------------------------------------------------------------------------
# Helper: write a NEXUS partition file.
# Usage: write_nexus <nexus_path> <model> <seqtype> <len1> <len2> ...
# -----------------------------------------------------------------------------
write_nexus() {
    local nex_path="$1"
    local model="$2"
    local seqtype="$3"
    shift 3
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
        local cp=""
        for ((i = 0; i < n; i++)); do
            [[ $i -gt 0 ]] && cp+=", "
            cp+="${model}:part$((i + 1)){1.0}"
        done
        echo "  charpartition myparts = $cp;"
        echo "end;"
    } > "$nex_path"
}

# -----------------------------------------------------------------------------
# Helper: run AliSim for one scenario.
# Usage: run_alisim <scenario_name> <model> <seqtype> <len1> <len2> ...
# Appends one row to META_TSV.
# -----------------------------------------------------------------------------
run_alisim() {
    local sname="$1"
    local model="$2"
    local seqtype="$3"
    shift 3
    local lengths=("$@")
    local n_parts=${#lengths[@]}
    local total_len=0
    for l in "${lengths[@]}"; do (( total_len += l )); done
    local min_len=${lengths[0]} max_len=${lengths[0]}
    for l in "${lengths[@]}"; do
        (( l < min_len )) && min_len=$l
        (( l > max_len )) && max_len=$l
    done

    local prefix="$OUTDIR/$sname"
    local nex_out="${prefix}_partition.nex"
    local phy_out="${prefix}.phy"

    # Codon charset positions are in nucleotides — multiply each length by 3
    local nex_lengths=("${lengths[@]}")
    if [[ "$seqtype" == "CODON" ]]; then
        for ((i = 0; i < ${#nex_lengths[@]}; i++)); do
            nex_lengths[i]=$(( nex_lengths[i] * 3 ))
        done
    fi

    write_nexus "$nex_out" "$model" "$seqtype" "${nex_lengths[@]}"

    local err_tmp
    err_tmp=$(mktemp)

    "$IQTREE" --alisim "$prefix" \
        -p "$nex_out" \
        -t "RANDOM{yh/$NTAXA}" \
        -st "$seqtype" \
        --seed "$SEED" \
        --redo -quiet 2>"$err_tmp"
    local status=$?

    if [[ $status -eq 0 && -f "$phy_out" ]]; then
        printf "  OK   %-35s  %2d parts  total %6d sites\n" \
            "$sname" "$n_parts" "$total_len"
        printf "%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n" \
            "$sname" "$seqtype" "$n_parts" "$total_len" "$min_len" "$max_len" \
            "$phy_out" "$nex_out" >> "$META_TSV"
        rm -f "$err_tmp"
        return 0
    else
        printf "  FAIL %-35s\n" "$sname"
        if [[ -s "$err_tmp" ]]; then
            echo "    --- IQ-TREE error output ---"
            sed 's/^/    /' "$err_tmp"
            echo "    ---"
        fi
        rm -f "$err_tmp"
        printf "%s\t%s\t%d\t%d\t%d\t%d\tNA\tNA\n" \
            "$sname" "$seqtype" "$n_parts" "$total_len" "$min_len" "$max_len" >> "$META_TSV"
        return 1
    fi
}

# -----------------------------------------------------------------------------
# Write metadata TSV header
# -----------------------------------------------------------------------------
META_TSV="$OUTDIR/scenarios.tsv"
echo -e "scenario\tseq_type\tn_parts\ttotal_len\tmin_len\tmax_len\tphy\tnex" > "$META_TSV"

# -----------------------------------------------------------------------------
# Run all 6 scenarios for each data type
# -----------------------------------------------------------------------------
for dtype_idx in "${!DTYPE_NAMES[@]}"; do
    dtype="${DTYPE_NAMES[$dtype_idx]}"
    model="${DTYPE_MODELS[$dtype_idx]}"

    echo "=== $dtype  model=$model  very_light=${VERY_LIGHT}  very_heavy=${VERY_HEAVY} ==="

    # a: 2 very light
    run_alisim "${dtype}_a_2_very_light"   "$model" "$dtype" $VERY_LIGHT $VERY_LIGHT

    # b: 10 very light
    b_lens=(); for ((i=0; i<10; i++)); do b_lens+=($VERY_LIGHT); done
    run_alisim "${dtype}_b_10_very_light"  "$model" "$dtype" "${b_lens[@]}"

    # c: 2 very heavy
    run_alisim "${dtype}_c_2_very_heavy"   "$model" "$dtype" $VERY_HEAVY $VERY_HEAVY

    # d: 10 very heavy
    d_lens=(); for ((i=0; i<10; i++)); do d_lens+=($VERY_HEAVY); done
    run_alisim "${dtype}_d_10_very_heavy"  "$model" "$dtype" "${d_lens[@]}"

    # e: 2 mixed (1 very heavy + 1 very light)
    run_alisim "${dtype}_e_2_mixed"        "$model" "$dtype" $VERY_HEAVY $VERY_LIGHT

    # f: 10 mixed (alternating heavy/light: 10002 102 10002 102 ...)
    f_lens=()
    for ((i=0; i<5; i++)); do f_lens+=($VERY_HEAVY $VERY_LIGHT); done
    run_alisim "${dtype}_f_10_mixed"       "$model" "$dtype" "${f_lens[@]}"

    echo ""
done

echo "Metadata saved to: $META_TSV"
echo ""
echo "Scenario summary:"
column -t -s $'\t' "$META_TSV"
