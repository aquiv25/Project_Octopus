"""
G4 quadruplex pattern search in promoters across cephalopod species.

Pipeline:
1. Parse GFF → find target gene TSS coordinates (using product= descriptions in CDS)
2. Fetch promoter sequences via NCBI efetch (2kb upstream of TSS)
3. Search for G4 pattern (G3+N1-7 x4) on both strands
4. Parse G4 structure (G-runs, loops, bulges)
5. Save table + FASTA for alignment
"""
import re, gzip, json, time, subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

# ── Paths ────────────────────────────────────────────────────────────────
BASE    = Path("/Users/dassagaripova/Documents/Claude/Projects/Article/g4_promoter")
GFF_DIR = BASE / "gff"
OUT_DIR = BASE / "promoter_seqs"
OUT_DIR.mkdir(exist_ok=True)

PROMOTER_UP   = 2000   # bp upstream of TSS
PROMOTER_DOWN = 200    # bp downstream

# ── Target genes ─────────────────────────────────────────────────────────
KEYWORDS = {
    "DNMT1": ["dna methyltransferase 1", "dna (cytosine-5)-methyltransferase 1",
              "dna (cytosine-5) methyltransferase 1", "dnmt1"],
    "TET3":  ["methylcytosine dioxygenase tet3", "tet3", "tet methylcytosine dioxygenase"],
    "EZH2":  ["histone-lysine n-methyltransferase ezh2", "enhancer of zeste homolog 2", "ezh2"],
    "KMT2A": ["histone-lysine n-methyltransferase 2a", "lysine methyltransferase 2a",
              "mixed-lineage leukemia", "kmt2a"],
    "CHD1":  ["chromodomain-helicase-dna-binding protein 1", "chd1"],
    "BRD3":  ["bromodomain-containing protein 3", "brd3"],
    "ADAR1": ["double-stranded rna-specific adenosine deaminase",
              "adenosine deaminase acting on rna", "adar"],
    "UHRF1": ["e3 ubiquitin-protein ligase uhrf1", "uhrf1",
              "ubiquitin-like with phd and ring finger domains"],
}

# ── Species GFF files ─────────────────────────────────────────────────────
SPECIES_GFF = {
    "O_vulgaris": GFF_DIR / "O_vulgaris.gff.gz",
    "O_sinensis": GFF_DIR / "O_sinensis.gff.gz",
    "A_fangsiao": GFF_DIR / "A_fangsiao.gff.gz",
    "A_hians":    GFF_DIR / "A_hians.gff.gz",
}

# ── G4 patterns ──────────────────────────────────────────────────────────
# Classic parallel G4: G3+N1-7G3+N1-7G3+N1-7G3+
G4_RE = re.compile(r'(G{3,6}[ACGT]{1,7}G{3,6}[ACGT]{1,7}G{3,6}[ACGT]{1,7}G{3,6})', re.IGNORECASE)
# Also antiparallel on C-strand
G4_RC_RE = re.compile(r'(C{3,6}[ACGT]{1,7}C{3,6}[ACGT]{1,7}C{3,6}[ACGT]{1,7}C{3,6})', re.IGNORECASE)


# ════════════════════════════════════════════════════════════════════════
# Step 1: Parse GFF files
# ════════════════════════════════════════════════════════════════════════
def parse_gff(gff_path, sp_name):
    """Find best CDS record for each target gene. Returns dict gene->{chrom,tss,strand,desc}."""
    found = {}
    print(f"\n  Parsing {sp_name} GFF...")
    with gzip.open(gff_path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9: continue
            if parts[2] not in ("CDS", "mRNA", "gene"): continue
            chrom = parts[0]
            start, end = int(parts[3]), int(parts[4])
            strand = parts[6]
            attrs = parts[8]

            # Extract product/name
            desc = ""
            for key in ("product=", "Name=", "gene="):
                m = re.search(rf'{key}([^;]+)', attrs)
                if m:
                    desc += " " + m.group(1).lower()

            for gene, kws in KEYWORDS.items():
                if gene not in found:
                    if any(kw.lower() in desc for kw in kws):
                        # TSS: start if +, end if -
                        tss = start if strand == "+" else end
                        found[gene] = {
                            "chrom": chrom, "tss": tss,
                            "strand": strand, "gene_start": start,
                            "gene_end": end, "desc": desc.strip()
                        }
    return found


# ════════════════════════════════════════════════════════════════════════
# Step 2: Fetch promoter sequences via NCBI efetch
# ════════════════════════════════════════════════════════════════════════
def fetch_promoter(chrom, tss, strand, gene_name, sp_name):
    """Fetch 2kb upstream sequence via NCBI eutils REST API."""
    import urllib.request, urllib.parse

    if strand == "+":
        seq_start = max(1, tss - PROMOTER_UP)
        seq_stop  = tss + PROMOTER_DOWN
    else:
        seq_start = max(1, tss - PROMOTER_DOWN)
        seq_stop  = tss + PROMOTER_UP

    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={chrom}"
        f"&seq_start={seq_start}&seq_stop={seq_stop}"
        "&rettype=fasta&retmode=text&email=researcher@example.com"
    )
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            text = resp.read().decode("utf-8")
        if text.startswith(">"):
            lines = text.strip().split("\n")
            seq = "".join(lines[1:]).upper()
            if strand == "-":
                seq = str(Seq(seq).reverse_complement())
            return seq
    except Exception as e:
        pass
    return None


# ════════════════════════════════════════════════════════════════════════
# Step 3: G4 pattern search
# ════════════════════════════════════════════════════════════════════════
def find_g4s(seq, strand_label="+"):
    """Find all G4 motifs in sequence. Returns list of dicts."""
    results = []
    for m in G4_RE.finditer(seq):
        results.append({"seq": m.group(), "pos": m.start(),
                        "strand": strand_label, "type": "G4"})
    # Check C-strand (antiparallel)
    rc = str(Seq(seq).reverse_complement())
    for m in G4_RE.finditer(rc):
        results.append({"seq": m.group(), "pos": len(seq) - m.end(),
                        "strand": "-" if strand_label == "+" else "+",
                        "type": "G4_antisense"})
    return results


def parse_g4_structure(seq):
    """Parse G4 into G-runs and loops. Returns structural dict."""
    m = re.match(r'(G{2,6})([ACGT]{1,12})(G{2,6})([ACGT]{1,12})(G{2,6})([ACGT]{1,12})(G{2,6})',
                 seq, re.IGNORECASE)
    if not m:
        return None
    g1, l1, g2, l2, g3, l3, g4 = m.groups()
    return {
        "sequence":     seq,
        "G_runs":       [g1, g2, g3, g4],
        "G_run_lens":   [len(g1), len(g2), len(g3), len(g4)],
        "loops":        [l1, l2, l3],
        "loop_lengths": [len(l1), len(l2), len(l3)],
        "total_len":    len(seq),
        "has_bulge":    any(len(g) > 4 for g in [g1, g2, g3, g4]),
        "G_run_uniform": len(set([len(g1), len(g2), len(g3), len(g4)])) == 1,
        "loop_sum":     len(l1) + len(l2) + len(l3),
    }


# ════════════════════════════════════════════════════════════════════════
# Main pipeline
# ════════════════════════════════════════════════════════════════════════
all_results = {}   # sp_name → {gene → {coords, g4_list}}

print("=" * 65)
print("STEP 1: Find gene TSS from GFF files")
print("=" * 65)

gene_coords = {}
for sp_name, gff_path in SPECIES_GFF.items():
    if not gff_path.exists():
        print(f"  ⚠  {sp_name}: GFF not found")
        continue
    coords = parse_gff(gff_path, sp_name)
    gene_coords[sp_name] = coords
    print(f"  {sp_name}: {len(coords)}/8 genes → {sorted(coords.keys())}")

with open(BASE / "gene_coords.json", "w") as f:
    json.dump(gene_coords, f, indent=2)

print("\n" + "=" * 65)
print("STEP 2: Fetch promoter sequences via NCBI efetch")
print("=" * 65)

all_promoters = {}   # sp_name → {gene → seq}
for sp_name, coords in gene_coords.items():
    all_promoters[sp_name] = {}
    for gene, info in coords.items():
        print(f"  {sp_name} / {gene}: chr={info['chrom']}, TSS={info['tss']}, strand={info['strand']}")
        seq = fetch_promoter(info["chrom"], info["tss"], info["strand"], gene, sp_name)
        if seq:
            all_promoters[sp_name][gene] = seq
            print(f"    ✅ fetched {len(seq)} bp")
        else:
            print(f"    ❌ fetch failed")
        time.sleep(0.35)  # NCBI rate limit

# Save promoter FASTAs
prom_fasta = BASE / "promoters_all_species.fasta"
with open(prom_fasta, "w") as f:
    for sp, genes in all_promoters.items():
        for gene, seq in genes.items():
            f.write(f">{gene}_{sp}_promoter\n{seq}\n")
print(f"\nSaved {prom_fasta}")

print("\n" + "=" * 65)
print("STEP 3: G4 pattern search in promoters")
print("=" * 65)

all_g4 = []   # list of dicts for table

for sp_name, genes in all_promoters.items():
    for gene, seq in genes.items():
        hits = find_g4s(seq)
        for hit in hits:
            struct = parse_g4_structure(hit["seq"])
            if struct:
                row = {
                    "species": sp_name,
                    "gene": gene,
                    **struct,
                    "position_in_promoter": hit["pos"],
                    "dist_to_tss": 2000 - hit["pos"],
                }
                all_g4.append(row)

print(f"  Total G4 motifs found: {len(all_g4)}")
for sp in SPECIES_GFF:
    n = sum(1 for r in all_g4 if r["species"] == sp)
    genes_with = set(r["gene"] for r in all_g4 if r["species"] == sp)
    print(f"  {sp}: {n} G4s in {len(genes_with)} genes: {sorted(genes_with)}")

# Save table
import csv
table_path = BASE / "g4_in_promoters.tsv"
if all_g4:
    with open(table_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_g4[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(all_g4)
    print(f"\nSaved table: {table_path}")

# Save G4 sequences for alignment
g4_fasta = BASE / "g4_sequences.fasta"
with open(g4_fasta, "w") as f:
    for i, row in enumerate(all_g4):
        f.write(f">{row['gene']}_{row['species']}_g4_{i}\n{row['sequence']}\n")
print(f"Saved G4 FASTA: {g4_fasta}")

print("\n✅ Done!")
