# Example Data Formats

This directory documents the expected input file formats.

## enrichment.all.tsv
STRING enrichment export (TSV, tab-separated). Download from string-db.org:
`Your proteins → Analysis → Functional Enrichment → Export (All enrichment)`

Columns (no header, `#` comment lines):
```
category | term_id | description | obs | bg | strength | signal | fdr | ids | labels
```
Example:
```
GO Process	GO:0006355	regulation of DNA-templated transcription	42	180	1.23	2.1	0.0001	...	GENE1,GENE2,...
```

## genes_with_PQS_in_promoters.txt
One gene name per line (output of `g4_promoter_analysis.py`):
```
ARID1A
ARID1B
ARID2
...
```

## zdna_thr025.bedgraph
4-column BEDGraph (ZHunter output):
```
NC_024355.1   1000   1016   0.87
NC_024355.1   2400   2416   1.23
```
