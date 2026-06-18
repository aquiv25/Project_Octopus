#!/bin/bash
# MAFFT alignment + IQ-TREE phylogenetic trees for 8 epigenetic genes
set -e

BASE="/Users/dassagaripova/Documents/Claude/Projects/Article"
GENES_DIR="$BASE/per_gene_mega"
ALIGNED_DIR="$BASE/aligned"
TREES_DIR="$BASE/trees"

mkdir -p "$ALIGNED_DIR" "$TREES_DIR"

echo "============================================"
echo "STEP 1: MAFFT alignments"
echo "============================================"

for fasta in "$GENES_DIR"/*_8species.fasta; do
    gene=$(basename "$fasta" _8species.fasta)
    n=$(grep -c ">" "$fasta")
    echo ""
    echo "▶ $gene ($n sequences)"

    if [ "$n" -lt 2 ]; then
        echo "  ⚠ Skipping — need at least 2 sequences"
        continue
    fi

    out="$ALIGNED_DIR/${gene}_aligned.fasta"
    mafft --auto --reorder "$fasta" > "$out" 2>/dev/null
    n_aln=$(grep -c ">" "$out")
    echo "  ✅ Aligned: $n_aln sequences → $(basename $out)"
done

echo ""
echo "============================================"
echo "STEP 2: IQ-TREE phylogenetic trees"
echo "============================================"

for aln in "$ALIGNED_DIR"/*_aligned.fasta; do
    gene=$(basename "$aln" _aligned.fasta)
    n=$(grep -c ">" "$aln")

    if [ "$n" -lt 3 ]; then
        echo "⚠ $gene: need ≥3 sequences for tree (have $n), skipping"
        continue
    fi

    echo ""
    echo "▶ $gene ($n sequences)"
    prefix="$TREES_DIR/${gene}_tree"

    iqtree -s "$aln" -m LG+G4 -B 1000 -T AUTO \
           --prefix "$prefix" --redo -quiet 2>/dev/null

    if [ -f "${prefix}.treefile" ]; then
        echo "  ✅ Tree → ${gene}_tree.treefile"
        echo "  $(cat ${prefix}.treefile)"
    else
        echo "  ❌ Tree not generated"
    fi
done

echo ""
echo "============================================"
echo "All done! Files:"
echo "  Alignments: $ALIGNED_DIR"
echo "  Trees:      $TREES_DIR"
echo "============================================"
