#!/bin/bash
# =============================================================================
# simulate_alignments.sh
#
# Simulate alignments of varying sizes and data types using IQ-TREE AliSim,
# for benchmarking ModelFinder thread-capping behaviour.
#
# 8 alignment lengths (k sites/codons) × 3 data types = 24 datasets.
#
# k values: 51 102 201 501 1002 2001 5001 10002
# All are divisible by 3, so CODON --length (k*3 nucleotides) is always valid.
#
# Data types:
#   DNA   (GTR+G, 4  states)
#   AA    (LG+G,  20 states)
#   CODON (GY+G,  61 states)  --length passed as k*3 nucleotides
#
# Usage:
#   bash simulate_alignments.sh [IQTREE_BIN] [OUTDIR]
# =============================================================================

IQTREE=${1:-iqtree3}
OUTDIR=${2:-sim_alignments}
NTAXA=20
SEED=12345

mkdir -p "$OUTDIR"

if ! command -v "$IQTREE" &>/dev/null && [[ ! -x "$IQTREE" ]]; then
    echo "ERROR: IQ-TREE binary not found: $IQTREE"
    exit 1
fi

echo "IQ-TREE binary   : $IQTREE"
echo "Output directory : $OUTDIR"
echo "Taxa: $NTAXA  Seed: $SEED"
echo ""

# k = alignment length in sites (codons for CODON)
K_VALUES=(51 102 201 501 1002 2001 5001 10002)

TREE="RANDOM{yh/$NTAXA}"

run_alisim() {
    local prefix="$1"
    local model="$2"
    local length="$3"    # --length value passed to AliSim
    local seqtype="$4"

    "$IQTREE" --alisim "$prefix" \
        -m "$model" \
        -t "$TREE" \
        --length "$length" \
        --seed "$SEED" \
        -st "$seqtype" \
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
for k in "${K_VALUES[@]}"; do
    run_alisim "$OUTDIR/dna_n${NTAXA}_k${k}" "GTR+G" "$k" "DNA"
done
echo ""

# =============================================================================
# Protein alignments  (LG+G, 20 states)
# =============================================================================
echo "=== Protein alignments (LG+G, 20 states) ==="
for k in "${K_VALUES[@]}"; do
    run_alisim "$OUTDIR/aa_n${NTAXA}_k${k}" "LG+G" "$k" "AA"
done
echo ""

# =============================================================================
# Codon alignments  (GY+G, 61 states)
# --length is in nucleotides (k * 3); all k are divisible by 3.
# =============================================================================
echo "=== Codon alignments (GY+G, 61 states) ==="
for k in "${K_VALUES[@]}"; do
    nt_len=$(( k * 3 ))
    run_alisim "$OUTDIR/codon_n${NTAXA}_k${k}" "GY+G" "$nt_len" "CODON"
done
echo ""

echo "Simulation complete. Alignments saved to: $OUTDIR/"
echo ""
echo "File summary:"
ls "$OUTDIR"/*.phy 2>/dev/null | while read f; do
    sites=$(awk 'NR==1{print $2}' "$f" 2>/dev/null)
    echo "  $f  ($sites sites)"
done
