#!/bin/bash
# =============================================================================
# simulate_partitions.sh
#
# Simulate multi-partition alignments using IQ-TREE AliSim for benchmarking
# the new thread-scheduling strategy in PartitionFinder.
#
# Three data types are simulated (DNA, AA, CODON), each with five partition
# scenarios that exercise different regions of the thread-cap formula
# max(1, nptn * nstate / 4000):
#
#   DNA  (4  states): 500 sites -> m_p=1,  2000 -> m_p=2,  8000 -> m_p=8
#   AA   (20 states): 500 sites -> m_p=2,  2000 -> m_p=10, 8000 -> m_p=40
#   CODON(61 states): 100 sites -> m_p=1,  500  -> m_p=7,  2000 -> m_p=30
#
# Scenarios (site counts are in alignment sites: bp for DNA/AA, codons for CODON):
#   A – Many uniform light   (16 x 500)
#   B – Mixed heavy + light  (1 x 8000 + 15 x 200)
#   C – Few heavy            (4 x 4000)
#   D – Many heavy           (8 x 2000)
#   E – Very many very light (32 x 100)
#
# Usage:
#   bash simulate_partitions.sh [IQTREE_BIN] [OUTDIR]
# =============================================================================

IQTREE=${1:-iqtree3}
OUTDIR=${2:-sim_partitions}
NTAXA=20
SEED=12345

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
# Data types: parallel arrays (index-aligned)
# -----------------------------------------------------------------------------
DTYPE_NAMES=("DNA"   "AA"   "CODON")
DTYPE_MODELS=("GTR+G" "LG+G" "GY+G")

# -----------------------------------------------------------------------------
# Helper: write a NEXUS partition file.
# Usage: write_nexus <nexus_path> <model> <len1> <len2> ...
# The tree length {1.0} is required by AliSim under proportional branch lengths.
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
        # {1.0} = tree length per partition, required by AliSim (-p mode)
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
# Sets globals PHY_OUT and NEX_OUT. Appends one row to META_TSV.
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
    local min_len=${lengths[0]}
    local max_len=${lengths[0]}
    for l in "${lengths[@]}"; do
        (( l < min_len )) && min_len=$l
        (( l > max_len )) && max_len=$l
    done

    local prefix="$OUTDIR/$sname"
    NEX_OUT="${prefix}_partition.nex"
    PHY_OUT="${prefix}.phy"

    # Codon charset positions are in nucleotides — multiply each length by 3
    local nex_lengths=("${lengths[@]}")
    if [[ "$seqtype" == "CODON" ]]; then
        for ((i = 0; i < ${#nex_lengths[@]}; i++)); do
            nex_lengths[i]=$(( nex_lengths[i] * 3 ))
        done
    fi

    write_nexus "$NEX_OUT" "$model" "${nex_lengths[@]}"

    local err_tmp
    err_tmp=$(mktemp)

    "$IQTREE" --alisim "$prefix" \
        -p "$NEX_OUT" \
        -t "RANDOM{yh/$NTAXA}" \
        -st "$seqtype" \
        --seed "$SEED" \
        --redo -quiet 2>"$err_tmp"
    local status=$?

    if [[ $status -eq 0 && -f "$PHY_OUT" ]]; then
        printf "  OK   %-30s  %2d partitions  total %6d sites\n" \
            "$sname" "$n_parts" "$total_len"
        printf "%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n" \
            "$sname" "$seqtype" "$n_parts" "$total_len" "$min_len" "$max_len" \
            "$PHY_OUT" "$NEX_OUT" >> "$META_TSV"
        rm -f "$err_tmp"
        return 0
    else
        printf "  FAIL %-30s\n" "$sname"
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
# Build scenario length arrays (sites: bp for DNA/AA, codons for CODON)
# -----------------------------------------------------------------------------
A_LENS=();  for ((i=0; i<16; i++)); do A_LENS+=(500);  done   # 16 x 500
B_LENS=(8000); for ((i=0; i<15; i++)); do B_LENS+=(200); done # 1 x 8000 + 15 x 200
C_LENS=();  for ((i=0; i<4;  i++)); do C_LENS+=(4000); done   # 4 x 4000
D_LENS=();  for ((i=0; i<8;  i++)); do D_LENS+=(2000); done   # 8 x 2000
E_LENS=();  for ((i=0; i<32; i++)); do E_LENS+=(100);  done   # 32 x 100

# -----------------------------------------------------------------------------
# Run all scenarios for each data type
# -----------------------------------------------------------------------------
for dtype_idx in "${!DTYPE_NAMES[@]}"; do
    dtype="${DTYPE_NAMES[$dtype_idx]}"
    model="${DTYPE_MODELS[$dtype_idx]}"

    echo "=== $dtype scenarios (model: $model) ==="
    run_alisim "${dtype}_A_uniform_light"   "$model" "$dtype" "${A_LENS[@]}"
    run_alisim "${dtype}_B_mixed"           "$model" "$dtype" "${B_LENS[@]}"
    run_alisim "${dtype}_C_few_heavy"       "$model" "$dtype" "${C_LENS[@]}"
    run_alisim "${dtype}_D_many_heavy"      "$model" "$dtype" "${D_LENS[@]}"
    run_alisim "${dtype}_E_many_very_light" "$model" "$dtype" "${E_LENS[@]}"
    echo ""
done

echo "Metadata saved to: $META_TSV"
echo ""
echo "Scenario summary:"
column -t -s $'\t' "$META_TSV"
