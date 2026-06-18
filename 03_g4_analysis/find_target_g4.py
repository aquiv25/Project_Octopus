"""
Find strict G4 in promoters of 5 target genes across all 5 species.
Genes: tropomyosin, MED28, SIK1, homeobox abdominal-A, zinc finger protein 271
"""
import re, gzip, time, urllib.request, json, csv
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.Seq import Seq

BASE = Path("/Users/dassagaripova/Documents/Claude/Projects/Article/g4_promoter")
G4_STRICT = re.compile(r'G{3,6}[ACGT]{1,7}G{3,6}[ACGT]{1,7}G{3,6}[ACGT]{1,7}G{3,6}', re.IGNORECASE)
ALL_SPECIES = ["O_bimaculoides","O_vulgaris","O_sinensis","A_fangsiao","A_hians"]

# Keywords for each target gene
TARGETS = {
    "TPM":    ["tropomyosin"],
    "MED28":  ["mediator of rna polymerase ii transcription subunit 28",
               "mediator complex subunit 28", "med28"],
    "SIK1":   ["salt-inducible kinase 1", "serine/threonine-protein kinase sik1",
               "sik1", "salt inducible kinase"],
    "ABDA":   ["homeobox protein abdominal-a", "homeobox abdominal",
               "abdominal-a", "abd-a"],
    "ZNF271": ["zinc finger protein 271", "znf271"],
}

def build_index(gff_path):
    idx = defaultdict(list)
    seen = set()
    with gzip.open(gff_path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#"): continue
            p = line.strip().split("\t")
            if len(p) < 9 or p[2] != "CDS": continue
            prod_m = re.search(r'product=([^;]+)', p[8])
            if not prod_m: continue
            prod = prod_m.group(1).strip().lower()
            chrom = p[0]; start, end = int(p[3]), int(p[4]); strand = p[6]
            tss = start if strand == "+" else end
            key = (prod, chrom, tss)
            if key not in seen:
                seen.add(key)
                idx[prod].append((chrom, tss, strand))
    return idx

print("Building GFF indices...", flush=True)
sp_idx = {sp: build_index(BASE / f"gff/{sp}.gff.gz") for sp in ALL_SPECIES}
for sp in ALL_SPECIES:
    print(f"  {sp}: {len(sp_idx[sp])} products", flush=True)

def find_gene(gene_name, kws, idx):
    """Find best matching product in index."""
    for prod, coords in idx.items():
        for kw in kws:
            if kw in prod:
                return prod, coords[0]
    return None, None

def fetch_promoter(chrom, tss, strand, up=3000, down=300):
    if strand == "+": s = max(1, tss-up); e = tss+down
    else: s = max(1, tss-down); e = tss+up
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={chrom}&seq_start={s}&seq_stop={e}"
           f"&rettype=fasta&retmode=text&email=researcher@example.com")
    try:
        with urllib.request.urlopen(url, timeout=15) as r:
            text = r.read().decode()
        if not text.startswith(">"): return None
        seq = "".join(text.strip().split("\n")[1:]).upper()
        if strand == "-": seq = str(Seq(seq).reverse_complement())
        return seq
    except Exception as e:
        print(f"    fetch error: {e}", flush=True)
        return None

def best_g4(seq):
    hits = list(G4_STRICT.finditer(seq))
    rc = str(Seq(seq).reverse_complement())
    hits += list(G4_STRICT.finditer(rc))
    if not hits: return None
    # Return longest hit
    return max(hits, key=lambda m: len(m.group())).group()

# ── Main search ──────────────────────────────────────────────────────────────
print("\nSearching for G4 in target gene promoters...", flush=True)
rows = []

for gene, kws in TARGETS.items():
    print(f"\n=== {gene} ===", flush=True)
    for sp in ALL_SPECIES:
        prod_name, coords = find_gene(gene, kws, sp_idx[sp])
        if coords is None:
            print(f"  {sp}: NOT FOUND in GFF", flush=True)
            continue
        chrom, tss, strand = coords
        print(f"  {sp}: {prod_name[:50]} | {chrom} tss={tss} {strand}", flush=True)
        seq = fetch_promoter(chrom, tss, strand)
        time.sleep(0.35)
        if seq is None:
            print(f"    -> fetch failed", flush=True)
            continue
        g4 = best_g4(seq)
        if g4:
            print(f"    -> G4 found: {g4}", flush=True)
            rows.append({"gene": gene, "species": sp, "G4_seq": g4,
                         "type": "strict", "product": prod_name})
        else:
            print(f"    -> no strict G4 in promoter", flush=True)

print(f"\n{'='*50}", flush=True)
print(f"Total hits: {len(rows)}", flush=True)
for r in rows:
    print(f"  {r['gene']} / {r['species']}: {r['G4_seq']}", flush=True)

# Save
with open(BASE/"results/target_g4_strict.tsv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene","species","G4_seq","type","product"], delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"\nSaved to results/target_g4_strict.tsv", flush=True)
