# 07 — A-to-I Editing × Secondary DNA Structure Enrichment

Permutation-based enrichment analysis of Z-DNA and G-quadruplex (G4) structures
at A-to-I (ADAR) RNA editing sites in *Octopus bimaculoides*.

## Input data

| File | Source | Description |
|------|--------|-------------|
| `data/GSE284564_Oct_editing.csv` | GEO: GSE284564 | A-to-I editing sites (cold vs warm) |
| `02_zdna_analysis/results/zhunter_z-dna.bedgraph.gz` | ZHunter | Z-DNA predictions (thermodynamic) |
| `02_zdna_analysis/results/GCF_001194135.2_ASM119413v2_zdna_thr025.bedgraph.gz` | Z-DNABERT | Z-DNA predictions (deep learning, thr=0.25) |
| `03_g4_analysis/results/g4_pattern.bed` | Pattern search | G4/PQS predictions |

## Scripts

| Script | Description |
|--------|-------------|
| `scripts/zdna_ati_analysis.py` | Z-DNA (ZHunter) enrichment at editing sites |
| `scripts/zdna_ati_analysis_zdnabert.py` | Z-DNA (Z-DNABERT) enrichment at editing sites |
| `scripts/g4_ati_analysis.py` | G4 enrichment at editing sites |

## Method

For each set of editing sites (All A→G, Cold-regulated↑, Warm-regulated↑):
1. Count sites overlapping secondary structure regions within windows ±0, ±50, ±100, ±200, ±500 bp
2. Generate 1000 random permutations (same chromosomal proportions)
3. Compute fold enrichment (FE) and empirical p-value

## Results

| Structure | Window | FE (All) | FE (Cold↑) | FE (Warm↑) | p-value |
|-----------|--------|----------|------------|------------|---------|
| ZHunter Z-DNA | ±50bp | 1.21 | 1.17 | **1.29** | <0.001 |
| ZHunter Z-DNA | ±100bp | 1.19 | 1.18 | 1.20 | <0.001 |
| ZHunter Z-DNA | ±200bp | 1.13 | 1.13 | 1.11 | <0.001 |
| Z-DNABERT | ±50–200bp | 0.20–0.37 | — | — | depleted |
| G4 (pattern) | ±100bp | **1.25** | 1.24 | 1.24 | 0.001–0.035 |
| G4 (pattern) | ±200bp | 1.15 | 1.13 | 1.19 | 0.003–0.052 |

## Output files

| File | Description |
|------|-------------|
| `results/zdna_ati_results.tsv` | ZHunter enrichment table (n=1000 perm) |
| `results/zdnabert_ati_results.tsv` | Z-DNABERT enrichment table (n=100 perm) |
| `results/g4_ati_results.tsv` | G4 enrichment table (n=1000 perm) |
| `results/zdna_ati_enrichment.png` | ZHunter enrichment figure |
| `results/zdnabert_ati_enrichment.png` | Z-DNABERT enrichment figure |
| `results/g4_ati_enrichment.png` | G4 enrichment figure |
