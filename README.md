# Project Octopus: Secondary DNA Structures in *Octopus bimaculoides*

A comprehensive bioinformatic analysis of Z-DNA and G-quadruplex (G4) structures in the *Octopus bimaculoides* genome, including promoter annotation, ortholog analysis, functional enrichment, and phylogenetic inference for epigenetic regulatory genes.

---

## Repository Structure

```
Project_Octopus/
├── 01_genome_setup/        # Download genome, fetch gene sequences across species
├── 02_zdna_analysis/       # Z-DNA prediction (ZHunter + Z-DNABERT), ortholog search, MEME
├── 03_g4_analysis/         # G4/PQS prediction, promoter analysis, G4-DNABERT
├── 04_phylogenetics/       # MAFFT alignment + IQ-TREE phylogenetic trees
├── 05_enrichment/          # STRING enrichment bubble plots and heatmaps
├── 06_visualization/       # Ideogram, pie charts, Venn diagrams, figures
└── data/example/           # Example input file formats
```

---

## Analysis Overview

### 1. Genome Setup (`01_genome_setup/`)
- Download *O. bimaculoides* genome (GCF_001194135.2) from NCBI
- Fetch protein/CDS sequences for 8 epigenetic genes across 5–8 species via NCBI Entrez

### 2. Z-DNA Analysis (`02_zdna_analysis/`)
Two independent prediction methods were used and their overlap analysed:

- **ZHunter** — thermodynamic scoring; input BedGraph (`GCF_001194135.2_ASM119413v2_zdna_thr025.bedgraph`, threshold 0.25); identified **1,781,183** Z-DNA regions (1.91% of genome, mean length 23.9 bp)
- **Z-DNABERT** — deep learning model; identified **1,588,238** regions (10.48% of genome, mean length 147.6 bp); intersection with ZHunter: **681,135** regions

Downstream analysis (both methods):
- Find Z-DNA regions overlapping promoters (−2000 to +200 bp from TSS) using bedtools
- Extract orthologous Z-DNA promoter sequences across *Octopus* species
- Visualize MEME motif 1 (E = 8.6 × 10⁻²³) as sequence logo

**Key result:** 9,313 genes with Z-DNA in promoters; enriched for transcriptional regulation (GO) and signaling pathways (KEGG)

### 3. G4 / PQS Analysis (`03_g4_analysis/`)
- Pattern-based G4 search: `G{3,}N{1,7}G{3,}N{1,7}G{3,}N{1,7}G{3,}`
- Deep learning approach: G4-DNABERT fine-tuned model on octopus genome
- Promoter overlap → 626 genes with G4 in promoters
- Ortholog G4 search for target genes (SIK1, ZNF271, MED28) across 5 species
- MEME motif analysis (E = 8.4 × 10⁻³)

**Key result:** G4 enriched 1.97× in promoters; 1.65× in 5′ UTR; depleted in CDS

### 4. Phylogenetics (`04_phylogenetics/`)
- MAFFT alignment (--auto --reorder) of 8 epigenetic gene families
- IQ-TREE tree inference: model LG+G4, 1000 ultrafast bootstraps
- Species: *O. bimaculoides*, *O. vulgaris*, *O. sinensis*, *A. fangsiao*, *A. hians*, *E. cirrhosa*

```bash
# Run full pipeline
bash run_align_trees.sh
```

### 5. Functional Enrichment (`05_enrichment/`)
- Input: STRING enrichment TSV (exported from string-db.org → Analysis → Enrichment)
- Bubble plot: −log10(FDR) × gene count, colored by category
- Heatmap: gene × GO-term matrix

### 6. Visualization (`06_visualization/`)
- Chromosomal ideogram (G4 / Z-DNA distribution across 30 chromosomes)
- Pie charts: genomic distribution of G4 and Z-DNA
- Venn diagrams: overlap between ZHunter and Z-DNABERT predictions; overlap between G4 pattern search and G4-DNABERT
- Octopus figure (figure panel assembly)

---

## Requirements

### Python 3.13.2

Install all dependencies:
```bash
pip install -r requirements.txt
```

| Package | Version |
|---------|---------|
| biopython | 1.87 |
| pandas | 3.0.3 |
| numpy | 2.4.6 |
| matplotlib | 3.10.9 |
| requests | 2.32.4 |
| scikit-learn | 1.9.0 |
| torch | 2.12.0 |
| transformers | 5.9.0 |

### R 4.5.2

| Package | Version |
|---------|---------|
| ggplot2 | 4.0.3 |
| ggseqlogo | 0.2.2 |

### External tools

| Tool | Version | Purpose |
|------|---------|---------|
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | 7.526 | Multiple sequence alignment (`--auto --reorder`) |
| [IQ-TREE](http://www.iqtree.org/) | 3.1.2 | Phylogenetic inference (model LG+G4, 1000 UFBoot) |
| [MEME Suite](https://meme-suite.org/) | 5.5.9 | De novo motif discovery (ZOOPS, width 6–20 nt) |
| [ZHunter](https://github.com/ikostits/ZHunter) | — | Thermodynamic Z-DNA prediction based on dinucleotide propensity scoring |
| [Z-DNABERT](https://github.com/Tsinghua-gongjing/Z-DNABERT) | — | BERT-based Z-DNA prediction; applied to *O. bimaculoides* genome (threshold 0.25) |
| [G4-DNABERT](https://github.com/Tsinghua-gongjing/DNABERT) | — | BERT-based model fine-tuned on G4-seq / G4-ChIP data; applied to octopus promoter windows |

---

## Data

### Included in this repository (via Git LFS)

| File | Location | Size |
|------|----------|------|
| `GCF_001194135.2_ASM119413v2_zdna_thr025.bedgraph.gz` | `02_zdna_analysis/results/` | 24 MB |
| `zhunter_z-dna.bedgraph.gz` | `02_zdna_analysis/results/` | 14 MB |
| `ZDNA_promoters-2.csv` | `02_zdna_analysis/results/` | 3.3 MB |
| `g4_pattern.bed` | `03_g4_analysis/results/` | 1.1 MB |
| `genes_with_PQS_in_promoters.txt` | `03_g4_analysis/results/` | 12 KB |
| `pqs_genes_obimaculoides.txt` | `03_g4_analysis/results/` | 28 KB |
| `pqs_in_promoters.tsv` | `03_g4_analysis/results/` | 56 KB |
| `enrichment_g4.tsv` | `05_enrichment/data/` | 52 KB |
| `enrichment_zdna.tsv` | `05_enrichment/data/` | 1.8 MB |

### Not included (too large — download separately)

| File | Source | Size | Notes |
|------|--------|------|-------|
| `GCF_001194135.2_ASM119413v2_genomic.fna` | NCBI RefSeq (GCF_001194135.2) | 2.2 GB | Reference genome assembly |
| `GCF_001194135.2_ASM119413v2_genomic.gff` | NCBI RefSeq (GCF_001194135.2) | ~50 MB | Genome annotation |
| Z-DNABERT BedGraph predictions | Z-DNABERT model output | ~500 MB | Whole-genome Z-DNA probability scores |
| G4-DNABERT model weights (`pytorch_model.bin`) | Fine-tuned DNABERT | 2 × ~360 MB | Stored in `03_g4_analysis/G4-DNABERT_model/` — download separately |

Download genome and annotation:
```bash
python 01_genome_setup/download_octopus_genome.py
```

---

## Citation

> Garipova D. et al. (2025). *Comprehensive bioinformatic analysis of Z-DNA and G-quadruplexes in the Octopus bimaculoides genome.* [Manuscript in preparation]

---

## Author

Daria Garipova — [GitHub](https://github.com/aquiv25)
