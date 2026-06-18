#!/bin/bash
set -e

BASE="/Users/dassagaripova/Documents/Claude/Projects/Article"
GENES_DIR="$BASE/per_gene_new"
ALIGNED_DIR="$BASE/aligned_new"
TREES_DIR="$BASE/trees_new"

mkdir -p "$ALIGNED_DIR" "$TREES_DIR"

echo "============================================"
echo "STEP 1: MAFFT alignments (new genes)"
echo "============================================"

for fasta in "$GENES_DIR"/*_species.fasta; do
    gene=$(basename "$fasta" _species.fasta)
    n=$(grep -c ">" "$fasta")
    echo ""
    echo "▶ $gene ($n sequences)"

    if [ "$n" -lt 2 ]; then
        echo "  ⚠ Skipping — need at least 2 sequences"
        continue
    fi

    out="$ALIGNED_DIR/${gene}_aligned.fasta"
    mafft --auto --reorder "$fasta" > "$out" 2>/dev/null
    echo "  ✅ Aligned → $(basename $out)"
done

echo ""
echo "============================================"
echo "STEP 2: IQ-TREE phylogenetic trees"
echo "============================================"

for aln in "$ALIGNED_DIR"/*_aligned.fasta; do
    gene=$(basename "$aln" _aligned.fasta)
    n=$(grep -c ">" "$aln")

    if [ "$n" -lt 3 ]; then
        echo "⚠ $gene: need ≥3 sequences (have $n), skipping"
        continue
    fi

    echo ""
    echo "▶ $gene ($n sequences)"
    prefix="$TREES_DIR/${gene}_tree"

    iqtree -s "$aln" -m LG+G4 -B 1000 -T AUTO \
           --prefix "$prefix" --redo -quiet 2>/dev/null

    if [ -f "${prefix}.treefile" ]; then
        echo "  ✅ Tree → ${gene}_tree.treefile"
    else
        echo "  ❌ Tree not generated"
    fi
done

echo ""
echo "Done! Trees in: $TREES_DIR"
