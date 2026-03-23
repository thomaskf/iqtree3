#!/bin/bash
# =============================================================================
# simulate_alignments.sh
#
# Simulate alignments of varying sizes and data types using IQ-TREE AliSim,
# for benchmarking ModelFinder thread-capping behaviour.
#
# Data types and their effective states:
#   DNA   : 4  states  -> cap triggers at nptn * 4  / 4000
#   AA    : 20 states  -> cap triggers at nptn * 20 / 4000
#   CODON : 60 states  -> cap triggers at nptn * 60 / 4000
#
# Usage:
#   bash simulate_alignments.sh [IQTREE_BIN] [OUTDIR]
#
# Examples:
#   bash simulate_alignments.sh
#   bash simulate_alignments.sh /usr/local/bin/iqtree3 ./sim_data
# =============================================================================

IQTREE=${1:-iqtree3}
OUTDIR=${2:-sim_alignments}
NTAXA=20          # number of taxa in the simulated tree
SEED=12345        # random seed for reproducibility
NREP=1            # number of replicates per setting

mkdir -p "$OUTDIR"

# Verify IQ-TREE binary
if ! command -v "$IQTREE" &>/dev/null; then
    echo "ERROR: IQ-TREE binary not found: $IQTREE"
    exit 1
fi

echo "IQ-TREE binary : $IQTREE"
echo "Output directory: $OUTDIR"
echo "Taxa: $NTAXA  Seed: $SEED"
echo ""

# -----------------------------------------------------------------------------
# Alignment lengths (number of sites) to simulate.
# Chosen so that nptn * nstate spans the range around the threshold formula
# nptn * nstate / 4000, giving effective thread counts of 1 through ~16.
#
#  DNA (4 states):   threshold sites ~ 1000 per thread
#  AA  (20 states):  threshold sites ~  200 per thread
#  Codon (60 states):threshold sites ~   67 per thread
# -----------------------------------------------------------------------------
DNA_LENGTHS=(200 500 1000 2000 4000 8000 16000)
AA_LENGTHS=(50  100  200  500  1000 2000  5000)
CODON_LENGTHS=(50  100  200  500  1000 2000  4000)  # in codons (x3 nucleotides)

# Tree topology shared across data types (Yule-Harding random tree)
TREE="RANDOM{yh/$NTAXA}"

run_alisim() {
    local prefix="$1"
    local model="$2"
    local length="$3"
    local seqtype="$4"   # DNA, AA, or CODON (empty = auto)

    local st_flag=""
    if [[ -n "$seqtype" ]]; then
        st_flag="-st $seqtype"
    fi

    "$IQTREE" --alisim "$prefix" \
        -m "$model" \
        -t "$TREE" \
        --length "$length" \
        --seed "$SEED" \
        $st_flag \
        --redo -quiet 2>/dev/null

    if [[ $? -eq 0 ]]; then
        echo "  OK  $prefix  (length=$length)"
    else
        echo "  FAIL $prefix"
    fi
}

# =============================================================================
# DNA alignments  (GTR+G, 4 states)
# =============================================================================
echo "=== DNA alignments (GTR+G, 4 states) ==="
for len in "${DNA_LENGTHS[@]}"; do
    prefix="$OUTDIR/dna_n${NTAXA}_l${len}"
    run_alisim "$prefix" "GTR+G" "$len" "DNA"
done
echo ""

# =============================================================================
# Protein alignments  (LG+G, 20 states)
# =============================================================================
echo "=== Protein alignments (LG+G, 20 states) ==="
for len in "${AA_LENGTHS[@]}"; do
    prefix="$OUTDIR/aa_n${NTAXA}_l${len}"
    run_alisim "$prefix" "LG+G" "$len" "AA"
done
echo ""

# =============================================================================
# Codon alignments  (GY+G, 60 states)
# --length for codon models specifies number of codons
# =============================================================================
echo "=== Codon alignments (GY+G, 60 states) ==="
for len in "${CODON_LENGTHS[@]}"; do
    prefix="$OUTDIR/codon_n${NTAXA}_l${len}"
    run_alisim "$prefix" "GY+G" "$len" "CODON"
done
echo ""

echo "Simulation complete. Alignments saved to: $OUTDIR/"
echo ""
echo "File summary:"
ls "$OUTDIR"/*.phy 2>/dev/null | while read f; do
    sites=$(awk 'NR==1{print $2}' "$f" 2>/dev/null)
    echo "  $f  ($sites sites)"
done
