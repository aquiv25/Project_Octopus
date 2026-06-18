"""
Find 5 orthologous genes with Z-DNA in promoters across:
  - O. bimaculoides (genes_with_ZDNA_in_promoters.txt -> LOC IDs)
  - O. vulgaris (selected_vulgaris_proteins.fasta -> protein names)
  - A. fangsiao (computed: search (CG)n/(CA·TG)n in promoters)

Z-DNA motif: alternating purine-pyrimidine, score by ZH-score proxy
"""
import re, gzip, time, urllib.request, csv
from pathlib import Path
from collections import defaultdict
from Bio.Seq import Seq

BASE  = Path("/Users/dassagaripova/Documents/Claude/Projects/Article/zdna_promoter")
BASE.mkdir(exist_ok=True)
GFF_DIR = Path("/Users/dassagaripova/Documents/Claude/Projects/Article/g4_promoter/gff")

# Z-DNA motif: runs of alternating purine-pyrimidine (CG, CA/TG, GC)
# Minimum 6 nt (3 dinucleotide repeats)
ZDNA_RE = re.compile(
    r'(?:(?:CG){3,}|(?:GC){3,}|(?:CA){3,}|(?:TG){3,}|(?:AT){3,}|(?:TA){3,}'
    r'|(?:[AG][CT]){4,}|(?:[CT][AG]){4,})', re.IGNORECASE)

# ── Step 1: O. bimaculoides — LOC IDs → product names ────────────────────────
print("=== O. bimaculoides: mapping LOC IDs to product names ===")
obi_locs = set(open('/var/folders/1l/9_k40zy95y7_ymxby1kwr6v40000gn/T/genes_with_ZDNA_in_promoters.txt').read().split())
print(f"  LOC IDs with Z-DNA: {len(obi_locs)}")

obi_loc2prod = {}
obi_loc2coords = {}
with gzip.open(GFF_DIR / "O_bimaculoides.gff.gz", "rt", errors="replace") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        if len(p) < 9 or p[2] != "CDS": continue
        attrs = p[8]
        gene_m = re.search(r'gene=(LOC\d+)', attrs)
        prod_m = re.search(r'product=([^;]+)', attrs)
        if gene_m and prod_m:
            loc = gene_m.group(1)
            if loc in obi_locs:
                prod = re.sub(r'%2C', ',', prod_m.group(1)).strip().lower()
                if loc not in obi_loc2prod and "uncharacterized" not in prod:
                    obi_loc2prod[loc] = prod
                    chrom = p[0]; start, end = int(p[3]), int(p[4]); strand = p[6]
                    tss = start if strand == "+" else end
                    obi_loc2coords[loc] = (chrom, tss, strand)

print(f"  Mapped named products: {len(obi_loc2prod)}")
obi_products = set(obi_loc2prod.values())

# ── Step 2: O. vulgaris — protein names from FASTA ───────────────────────────
print("\n=== O. vulgaris: extracting product names ===")
vul_products = set()
vul_prod2acc = {}
with open('/var/folders/1l/9_k40zy95y7_ymxby1kwr6v40000gn/T/selected_vulgaris_proteins.fasta') as f:
    for line in f:
        if line.startswith(">"):
            # >CAI9727365.1 mothers against decapentaplegic interacting [Octopus vulgaris]
            m = re.match(r'>(\S+)\s+(.+?)\s*\[', line)
            if m:
                acc = m.group(1)
                prod = m.group(2).strip().lower()
                if "hypothetical" not in prod:
                    prod_clean = re.sub(r'\s+isoform\s+x\d+', '', prod)
                    prod_clean = re.sub(r'\s+isoform\s+\d+', '', prod_clean)
                    vul_products.add(prod_clean)
                    vul_prod2acc[prod_clean] = acc

print(f"  O. vulgaris products with Z-DNA: {len(vul_products)}")

# ── Step 3: Find product names common to bimaculoides & vulgaris ─────────────
print("\n=== Finding common products (bimac ∩ vulgaris) ===")

def normalize(s):
    s = re.sub(r'\s+isoform\s+x\d+', '', s.lower())
    s = re.sub(r'\s+isoform\s+\d+', '', s)
    s = re.sub(r'-like$', '', s).strip()
    return s

obi_norm = {normalize(v): k for k, v in obi_loc2prod.items()}
vul_norm = {normalize(k): k for k in vul_products}

common_norm = set(obi_norm.keys()) & set(vul_norm.keys())
print(f"  Exact normalized matches: {len(common_norm)}")

# Partial match for remainder
partial = []
for on, loc in obi_norm.items():
    for vn in vul_norm:
        if on in vn or vn in on:
            if (on, vn) not in [(c, c) for c in common_norm]:
                partial.append((on, vn, loc))

print(f"  Partial matches: {len(partial)}")
all_common = list(common_norm)[:50]  # take up to 50 candidates

print(f"  Using {len(all_common)} candidates for A. fangsiao search")

# ── Step 4: A. fangsiao — find orthologs + compute Z-DNA ────────────────────
print("\n=== A. fangsiao: building GFF index ===")
fang_idx = defaultdict(list)
seen = set()
with gzip.open(GFF_DIR / "A_fangsiao.gff.gz", "rt", errors="replace") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        if len(p) < 9 or p[2] != "CDS": continue
        prod_m = re.search(r'product=([^;]+)', p[8])
        if not prod_m: continue
        prod = normalize(re.sub(r'%2C', ',', prod_m.group(1)))
        chrom = p[0]; start, end = int(p[3]), int(p[4]); strand = p[6]
        tss = start if strand == "+" else end
        key = (prod, chrom, tss)
        if key not in seen:
            seen.add(key)
            fang_idx[prod].append((chrom, tss, strand))

print(f"  A. fangsiao products: {len(fang_idx)}")

def fetch_promoter(chrom, tss, strand, up=2000, down=200):
    if strand == "+": s = max(1, tss-up); e = tss+down
    else: s = max(1, tss-down); e = tss+up
    url = (f"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={chrom}&seq_start={s}&seq_stop={e}"
           f"&rettype=fasta&retmode=text&email=researcher@example.com")
    for _ in range(2):
        try:
            with urllib.request.urlopen(url, timeout=15) as r:
                text = r.read().decode()
            if not text.startswith(">"): return None
            seq = "".join(text.strip().split("\n")[1:]).upper()
            if strand == "-": seq = str(Seq(seq).reverse_complement())
            return seq
        except: time.sleep(1)
    return None

def best_zdna(seq):
    hits = list(ZDNA_RE.finditer(seq))
    rc = str(Seq(seq).reverse_complement())
    hits += list(ZDNA_RE.finditer(rc))
    if not hits: return None
    best = max(hits, key=lambda m: len(m.group()))
    return best.group() if len(best.group()) >= 8 else None

# ── Step 5: Find genes with Z-DNA in all 3 species ───────────────────────────
print("\n=== Scanning for Z-DNA in all 3 species ===")
results = []

for norm_name in all_common:
    loc = obi_norm[norm_name]
    orig_prod = obi_loc2prod[loc]

    # Find in A. fangsiao
    fang_match = None
    for fn in fang_idx:
        if norm_name in fn or fn in norm_name:
            fang_match = fn; break
    if fang_match is None: continue

    # Fetch A. fangsiao promoter and check Z-DNA
    chrom, tss, strand = fang_idx[fang_match][0]
    seq = fetch_promoter(chrom, tss, strand)
    time.sleep(0.34)
    if seq is None: continue
    zdna_fang = best_zdna(seq)
    if zdna_fang is None:
        print(f"  skip {orig_prod[:40]} — no Z-DNA in A.fangsiao"); continue

    # Fetch O. vulgaris promoter (find coords from GFF)
    vul_gff = GFF_DIR / "O_vulgaris.gff.gz"
    vul_coords = None
    with gzip.open(vul_gff, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#"): continue
            p = line.strip().split("\t")
            if len(p) < 9 or p[2] != "CDS": continue
            prod_m = re.search(r'product=([^;]+)', p[8])
            if not prod_m: continue
            pn = normalize(re.sub(r'%2C', ',', prod_m.group(1)))
            if norm_name in pn or pn in norm_name:
                chrom2 = p[0]; start2, end2 = int(p[3]), int(p[4]); strand2 = p[6]
                tss2 = start2 if strand2 == "+" else end2
                vul_coords = (chrom2, tss2, strand2)
                break

    if vul_coords is None:
        print(f"  skip {orig_prod[:40]} — not in vulgaris GFF"); continue

    seq_vul = fetch_promoter(*vul_coords)
    time.sleep(0.34)
    if seq_vul is None: continue
    zdna_vul = best_zdna(seq_vul)
    if zdna_vul is None:
        print(f"  skip {orig_prod[:40]} — no Z-DNA in O.vulgaris"); continue

    # Fetch O. bimaculoides promoter
    obi_coords = obi_loc2coords[loc]
    seq_obi = fetch_promoter(*obi_coords)
    time.sleep(0.34)
    if seq_obi is None: continue
    zdna_obi = best_zdna(seq_obi)
    if zdna_obi is None:
        print(f"  skip {orig_prod[:40]} — no Z-DNA in O.bimaculoides"); continue

    results.append({
        "gene": orig_prod,
        "O_bimaculoides": zdna_obi,
        "O_vulgaris":     zdna_vul,
        "A_fangsiao":     zdna_fang,
    })
    print(f"  HIT: {orig_prod[:50]}")
    if len(results) >= 5:
        print("Got 5 genes, done!"); break

print(f"\nFinal: {len(results)} genes with Z-DNA in all 3 species")
for r in results:
    print(f"  {r['gene'][:50]}")
    for sp in ["O_bimaculoides","O_vulgaris","A_fangsiao"]:
        print(f"    {sp}: {r[sp]}")

# Save TSV
import json
rows = []
for r in results:
    for sp in ["O_bimaculoides","O_vulgaris","A_fangsiao"]:
        rows.append({"gene": r["gene"], "species": sp, "zdna_seq": r[sp]})

with open(BASE / "zdna_orthologs.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene","species","zdna_seq"], delimiter="\t")
    w.writeheader(); w.writerows(rows)
with open(BASE / "zdna_orthologs.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"Saved zdna_orthologs.tsv")
