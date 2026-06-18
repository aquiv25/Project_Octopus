"""
Full G4-in-promoters comparative analysis across 5 cephalopod species.

Strategy:
  1. Scan ALL O. bimaculoides promoters → genes with G4
  2. For those genes, fetch promoter seqs from other 4 species via NCBI
  3. G4 pattern search on all species
  4. Build per-gene G4 alignment + structural annotation
  5. Save table + FASTA + figure
"""
import re, gzip, json, time, urllib.request, csv, os
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
import matplotlib
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

BASE    = Path("/Users/dassagaripova/Documents/Claude/Projects/Article/g4_promoter")
GFF_DIR = BASE / "gff"
OUT_DIR = BASE / "results"
OUT_DIR.mkdir(exist_ok=True)

PROMOTER_UP   = 3000
PROMOTER_DOWN = 300

# ── G4 patterns ──────────────────────────────────────────────────────────
# Strict: G3+ tetrad (canonical)
G4_STRICT = re.compile(r'(G{3,6}[ACGT]{1,7}G{3,6}[ACGT]{1,7}G{3,6}[ACGT]{1,7}G{3,6})', re.IGNORECASE)
# Relaxed: G2+ (many biologically relevant G4s)
G4_RELAX  = re.compile(r'(G{2,6}[ACGT]{1,12}G{2,6}[ACGT]{1,12}G{2,6}[ACGT]{1,12}G{2,6})', re.IGNORECASE)

# ── Target genes ─────────────────────────────────────────────────────────
EPIG_KEYWORDS = {
    "DNMT1": ["dna methyltransferase 1", "cytosine-5)-methyltransferase 1"],
    "TET3":  ["tet3", "methylcytosine dioxygenase tet"],
    "EZH2":  ["histone-lysine n-methyltransferase ezh2", "enhancer of zeste"],
    "KMT2A": ["lysine methyltransferase 2a", "mixed-lineage leukemia"],
    "CHD1":  ["chromodomain-helicase-dna-binding protein 1"],
    "BRD3":  ["bromodomain-containing protein 3"],
    "ADAR1": ["double-stranded rna-specific adenosine deaminase", "adenosine deaminase"],
    "UHRF1": ["e3 ubiquitin-protein ligase uhrf1", "uhrf1", "ubiquitin-like with phd"],
}

# ── Species GFF files ─────────────────────────────────────────────────────
SPECIES_GFF = {
    "O_bimaculoides": GFF_DIR / "O_bimaculoides.gff.gz",
    "O_vulgaris":     GFF_DIR / "O_vulgaris.gff.gz",
    "O_sinensis":     GFF_DIR / "O_sinensis.gff.gz",
    "A_fangsiao":     GFF_DIR / "A_fangsiao.gff.gz",
    "A_hians":        GFF_DIR / "A_hians.gff.gz",
}

SP_LABELS = {
    "O_bimaculoides": "O. bimaculoides",
    "O_vulgaris":     "O. vulgaris",
    "O_sinensis":     "O. sinensis",
    "A_fangsiao":     "A. fangsiao",
    "A_hians":        "A. hians",
}
SP_COLORS = {
    "O_bimaculoides": "#D6604D",
    "O_vulgaris":     "#4393C3",
    "O_sinensis":     "#2CA02C",
    "A_fangsiao":     "#FF7F0E",
    "A_hians":        "#8C564B",
}

# ════════════════════════════════════════════════════════════════════════
def find_g4(seq, strict=True):
    """Find all G4 motifs on both strands. Returns list of seqs."""
    pat = G4_STRICT if strict else G4_RELAX
    hits = pat.findall(seq)
    rc   = pat.findall(str(Seq(seq).reverse_complement()))
    return hits + rc

def parse_g4_structure(seq):
    """Decompose G4 into G-runs and loops."""
    m = re.match(r'(G{2,7})([ACGTN]{1,15})(G{2,7})([ACGTN]{1,15})(G{2,7})([ACGTN]{1,15})(G{2,7})',
                 seq, re.IGNORECASE)
    if not m: return None
    g1, l1, g2, l2, g3, l3, g4 = m.groups()
    return {
        "sequence":   seq,
        "G1": g1, "L1": l1, "G2": g2, "L2": l2,
        "G3": g3, "L3": l3, "G4": g4,
        "g_runs":     [len(g1), len(g2), len(g3), len(g4)],
        "loop_lens":  [len(l1), len(l2), len(l3)],
        "loop_seqs":  [l1, l2, l3],
        "total_len":  len(seq),
        "has_bulge":  any(len(g) > 4 for g in [g1,g2,g3,g4]),
        "uniform_G":  len(set([len(g1),len(g2),len(g3),len(g4)])) == 1,
        "loop_sum":   len(l1)+len(l2)+len(l3),
    }

def fetch_seq(chrom, seq_start, seq_stop, strand="+", retries=3):
    """Fetch sequence from NCBI eutils."""
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={chrom}&seq_start={seq_start}&seq_stop={seq_stop}"
           "&rettype=fasta&retmode=text&email=researcher@example.com")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                text = r.read().decode()
            if text.startswith(">"):
                seq = "".join(text.strip().split("\n")[1:]).upper()
                return str(Seq(seq).reverse_complement()) if strand == "-" else seq
        except Exception as e:
            time.sleep(1)
    return None

def parse_gff(gff_path):
    """Return {gene_name: {chrom, tss, strand}} using CDS product annotations."""
    found = {}
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    with opener(gff_path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#"): continue
            p = line.strip().split("\t")
            if len(p) < 9 or p[2] not in ("CDS","mRNA"): continue
            chrom, start, end, strand = p[0], int(p[3]), int(p[4]), p[6]
            attrs = p[8]
            desc_m = re.search(r'product=([^;]+)', attrs)
            if not desc_m: continue
            text = desc_m.group(1).lower()
            for gene, kws in EPIG_KEYWORDS.items():
                if gene not in found and any(kw in text for kw in kws):
                    tss = start if strand == "+" else end
                    found[gene] = {"chrom": chrom, "tss": tss, "strand": strand}
    return found

# ════════════════════════════════════════════════════════════════════════
# STEP 1: O. bimaculoides — scan all promoters, find epigenetic gene G4s
# ════════════════════════════════════════════════════════════════════════
print("=" * 65)
print("STEP 1: O. bimaculoides — G4 in promoters (relaxed pattern)")
print("=" * 65)

# Load GFF to map LOC IDs to gene names
obi_gff = GFF_DIR / "O_bimaculoides.gff.gz"
obi_loc = {}   # gene_name → LOC_id
with gzip.open(obi_gff, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        if len(p) < 9 or p[2] not in ("CDS","mRNA"): continue
        attrs = p[8]
        desc_m = re.search(r'product=([^;]+)', attrs)
        parent_m = re.search(r'Parent=gene-([^;]+)', attrs)
        if not desc_m or not parent_m: continue
        text = desc_m.group(1).lower()
        loc  = parent_m.group(1)
        for gene, kws in EPIG_KEYWORDS.items():
            if gene not in obi_loc and any(kw in text for kw in kws):
                obi_loc[gene] = loc

print(f"  Gene → LOC: {obi_loc}")

# Scan promoters.fna for G4 in epigenetic gene promoters
PROM_FNA = Path("/Users/dassagaripova/@github.com/aquiv25/promoters.fna")
obi_g4 = {}   # gene → list of G4 seqs
obi_all_g4_genes = []  # all genes with G4 (for any analysis)

for rec in SeqIO.parse(str(PROM_FNA), "fasta"):
    seq = str(rec.seq).upper()
    hits = find_g4(seq, strict=False)
    if hits:
        obi_all_g4_genes.append(rec.id)
        # Check if it's an epigenetic gene
        for gene, loc in obi_loc.items():
            if loc in rec.id:
                obi_g4[gene] = hits[0]
                print(f"  ✅ {gene} ({loc}): {hits[0]}")

print(f"\n  Epigenetic genes with G4 (relaxed): {sorted(obi_g4.keys())}")
print(f"  Total promoters with G4: {len(obi_all_g4_genes)}")

# ════════════════════════════════════════════════════════════════════════
# STEP 2: Fetch promoters for all 4 other species for epigenetic genes
# ════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("STEP 2: Fetch epigenetic gene promoters from other 4 species")
print("=" * 65)

all_promoter_seqs = {"O_bimaculoides": {}}
# Add O_bimaculoides from existing data
for rec in SeqIO.parse(str(PROM_FNA), "fasta"):
    seq = str(rec.seq).upper()
    for gene, loc in obi_loc.items():
        if loc in rec.id:
            all_promoter_seqs["O_bimaculoides"][gene] = seq

for sp_name, gff_path in list(SPECIES_GFF.items())[1:]:  # skip O_bimaculoides
    all_promoter_seqs[sp_name] = {}
    print(f"\n  {sp_name}:")
    coords = parse_gff(gff_path)
    for gene, info in coords.items():
        chrom, tss, strand = info["chrom"], info["tss"], info["strand"]
        s_start = max(1, tss - PROMOTER_UP)   if strand == "+" else max(1, tss - PROMOTER_DOWN)
        s_stop  = tss + PROMOTER_DOWN          if strand == "+" else tss + PROMOTER_UP
        seq = fetch_seq(chrom, s_start, s_stop, strand)
        if seq:
            all_promoter_seqs[sp_name][gene] = seq
            hits = find_g4(seq, strict=False)
            print(f"    {gene}: {len(seq)} bp, {'✅ G4: ' + hits[0][:30] if hits else '— no G4'}")
        else:
            print(f"    {gene}: ❌ fetch failed")
        time.sleep(0.35)

# ════════════════════════════════════════════════════════════════════════
# STEP 3: G4 search in all promoters
# ════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("STEP 3: G4 search — all species × all genes")
print("=" * 65)

table_rows = []
fasta_seqs = []

for sp_name, genes in all_promoter_seqs.items():
    for gene, seq in genes.items():
        strict_hits = find_g4(seq, strict=True)
        relax_hits  = find_g4(seq, strict=False)
        for hit in (strict_hits or relax_hits)[:2]:  # keep best 2
            struct = parse_g4_structure(hit)
            if struct:
                row = {
                    "species":    sp_name,
                    "gene":       gene,
                    "G4_seq":     hit,
                    "type":       "strict" if hit in strict_hits else "relaxed",
                    "G_run_lengths": "-".join(map(str, struct["g_runs"])),
                    "loop1":      struct["L1"],
                    "loop2":      struct["L2"],
                    "loop3":      struct["L3"],
                    "loop1_len":  len(struct["L1"]),
                    "loop2_len":  len(struct["L2"]),
                    "loop3_len":  len(struct["L3"]),
                    "total_len":  struct["total_len"],
                    "has_bulge":  struct["has_bulge"],
                    "uniform_G":  struct["uniform_G"],
                }
                table_rows.append(row)
                label = f"{gene}_{sp_name}_g4"
                fasta_seqs.append((label, hit))
                print(f"  {sp_name:18s} {gene:6s}  {hit[:35]:35s}  G={row['G_run_lengths']}  L={row['loop1_len']},{row['loop2_len']},{row['loop3_len']}  {row['type']}")

print(f"\n  Total: {len(table_rows)} G4 motifs across {len(set(r['species'] for r in table_rows))} species")

# Save table
table_path = OUT_DIR / "g4_promoters_table.tsv"
if table_rows:
    with open(table_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=table_rows[0].keys(), delimiter="\t")
        w.writeheader(); w.writerows(table_rows)
print(f"  Saved table: {table_path}")

# Save G4 FASTA
g4_fasta = OUT_DIR / "g4_promoters.fasta"
with open(g4_fasta, "w") as f:
    for label, seq in fasta_seqs:
        f.write(f">{label}\n{seq}\n")
print(f"  Saved FASTA: {g4_fasta}")

# ════════════════════════════════════════════════════════════════════════
# STEP 4: Per-gene structural comparison figure
# ════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("STEP 4: Draw structural alignment figure")
print("=" * 65)

def draw_g4_alignment(gene_rows, gene_name, ax):
    """Draw color-coded G4 structural alignment for one gene."""
    ax.set_facecolor("white")
    if not gene_rows:
        ax.text(0.5, 0.5, f"{gene_name}\nno G4", ha="center", va="center",
                transform=ax.transAxes, color="#aaa", fontsize=11)
        ax.axis("off"); return

    y_labels = []
    for row_idx, row in enumerate(gene_rows):
        y = len(gene_rows) - row_idx - 1
        sp = row["species"]
        color = SP_COLORS.get(sp, "#333")

        # Parse structure
        struct = parse_g4_structure(row["G4_seq"])
        if not struct: continue

        x = 0
        parts = [
            (struct["G1"], "#1a6e2e"),  # G-run: dark green
            (struct["L1"], "#F4A460"),  # loop: sandy brown
            (struct["G2"], "#1a6e2e"),
            (struct["L2"], "#F4A460"),
            (struct["G3"], "#1a6e2e"),
            (struct["L3"], "#F4A460"),
            (struct["G4"], "#1a6e2e"),
        ]
        for seq_part, bg_color in parts:
            w = len(seq_part)
            rect = mpatches.FancyBboxPatch(
                (x, y - 0.42), w, 0.84,
                boxstyle="round,pad=0.05",
                linewidth=0, facecolor=bg_color, alpha=0.85
            )
            ax.add_patch(rect)
            for ci, ch in enumerate(seq_part):
                ax.text(x + ci + 0.5, y, ch, ha="center", va="center",
                        fontsize=7, fontweight="bold", color="white",
                        fontfamily="monospace")
            x += w

        y_labels.append(SP_LABELS.get(sp, sp))

    # Axis setup
    max_len = max(len(r["G4_seq"]) for r in gene_rows)
    ax.set_xlim(-1, max_len + 1)
    ax.set_ylim(-0.8, len(gene_rows) - 0.2)
    ax.set_yticks(range(len(gene_rows)))
    ax.set_yticklabels(reversed(y_labels), fontsize=9, style="italic")
    ax.set_title(gene_name, fontsize=12, fontweight="bold", pad=4)
    ax.set_xlabel("Position in G4 motif", fontsize=8)
    ax.tick_params(left=True, bottom=False, labelbottom=False)
    for spine in ["top","right","bottom"]: ax.spines[spine].set_visible(False)

    # Legend patches
    legend_elems = [
        mpatches.Patch(fc="#1a6e2e", label="G-run", ec="none"),
        mpatches.Patch(fc="#F4A460", label="Loop",  ec="none"),
    ]
    ax.legend(handles=legend_elems, loc="lower right", fontsize=7,
              frameon=False, ncol=2)

# Group by gene
from collections import defaultdict
by_gene = defaultdict(list)
for row in table_rows:
    by_gene[row["gene"]].append(row)

genes_with_g4 = sorted(by_gene.keys())
print(f"  Genes with G4: {genes_with_g4}")

if not genes_with_g4:
    print("  No G4 found — cannot draw figure")
else:
    ncols = min(3, len(genes_with_g4))
    nrows = (len(genes_with_g4) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 7, nrows * 3.5),
                             facecolor="white")
    if nrows * ncols == 1:
        axes = [[axes]]
    elif nrows == 1:
        axes = [axes]
    elif ncols == 1:
        axes = [[a] for a in axes]

    for idx, gene in enumerate(genes_with_g4):
        ax = axes[idx // ncols][idx % ncols]
        draw_g4_alignment(by_gene[gene], gene, ax)

    # Hide empty axes
    for idx in range(len(genes_with_g4), nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    # Species legend at bottom
    handles = [mpatches.Patch(fc=SP_COLORS[sp], label=SP_LABELS[sp], ec="none")
               for sp in SP_LABELS if sp in set(r["species"] for r in table_rows)]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               fontsize=10, frameon=False, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(
        "G4 quadruplex motifs in promoters of epigenetic genes — Cephalopoda\n"
        "G-runs shown in green, loops in orange  |  Pattern search G{2-6}N{1-12} × 4",
        fontsize=12, fontweight="bold", y=1.01
    )
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    out_svg = OUT_DIR / "g4_promoters_alignment.svg"
    out_png = OUT_DIR / "g4_promoters_alignment.png"
    fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
    fig.savefig(out_png, format="png", bbox_inches="tight", dpi=150, facecolor="white")
    print(f"  ✅ {out_svg}")
    print(f"  ✅ {out_png}")
    plt.close()

print("\n✅ All done!")
