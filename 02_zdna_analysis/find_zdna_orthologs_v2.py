"""
Find 5 orthologous genes with Z-DNA in promoters across:
  - A. fangsiao: LOC IDs in txt -> product names via O. bimaculoides GFF
  - O. vulgaris: protein names from FASTA
  - O. bimaculoides: Z-DNA bedgraph intersection with GFF gene promoters
"""
import re, gzip, time, urllib.request, csv, json
from pathlib import Path
from collections import defaultdict
from Bio.Seq import Seq

BASE    = Path("/Users/dassagaripova/Documents/Claude/Projects/Article/zdna_promoter")
GFF_DIR = Path("/Users/dassagaripova/Documents/Claude/Projects/Article/g4_promoter/gff")
BEDGRAPH = Path("/Users/dassagaripova/Downloads/project/zhunter_z-dna.bedgraph")

ZDNA_RE = re.compile(
    r'(?:(?:CG){3,}|(?:GC){3,}|(?:CA){4,}|(?:TG){4,}'
    r'|(?:[AG][CT]){4,}|(?:[CT][AG]){4,})', re.IGNORECASE)

def normalize(s):
    s = re.sub(r'%2[Cc]', ',', s)
    s = re.sub(r'\s+isoform\s+x\d+', '', s.lower())
    s = re.sub(r'\s+isoform\s+\d+', '', s)
    s = re.sub(r'-like$', '', s).strip()
    return s

# ── Step 1: A. fangsiao LOC IDs -> product names via O.bimaculoides GFF ──────
print("=== A. fangsiao: map LOC IDs -> product names ===")
fang_locs = set(open('/var/folders/1l/9_k40zy95y7_ymxby1kwr6v40000gn/T/genes_with_ZDNA_in_promoters.txt').read().split())
print(f"  LOC IDs: {len(fang_locs)}")

loc2prod = {}
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
            if loc in fang_locs:
                prod = normalize(prod_m.group(1))
                if loc not in loc2prod and "uncharacterized" not in prod:
                    loc2prod[loc] = prod

fang_products = set(loc2prod.values())
print(f"  Named products for A. fangsiao: {len(fang_products)}")

# ── Step 2: O. vulgaris protein names ────────────────────────────────────────
print("\n=== O. vulgaris: extract product names ===")
vul_products = set()
with open('/var/folders/1l/9_k40zy95y7_ymxby1kwr6v40000gn/T/selected_vulgaris_proteins.fasta') as f:
    for line in f:
        if line.startswith(">"):
            m = re.match(r'>(\S+)\s+(.+?)\s*\[', line)
            if m:
                prod = normalize(m.group(2))
                if "hypothetical" not in prod and len(prod) > 5:
                    vul_products.add(prod)
print(f"  O. vulgaris products: {len(vul_products)}")

# ── Step 3: Common genes between A.fangsiao and O.vulgaris ───────────────────
print("\n=== Intersecting A.fangsiao ∩ O.vulgaris ===")
common = fang_products & vul_products
print(f"  Exact matches: {len(common)}")

# Partial matches for more candidates
extra = []
for fp in fang_products:
    for vp in vul_products:
        if fp != vp and (fp in vp or vp in fp) and len(fp) > 10:
            extra.append(fp); break
all_candidates = list(common) + [e for e in extra if e not in common]
print(f"  Total candidates (exact + partial): {len(all_candidates)}")

# ── Step 4: O. bimaculoides Z-DNA genes from bedgraph ────────────────────────
print("\n=== O. bimaculoides: finding genes with Z-DNA in promoters ===")
# Load bedgraph regions (score threshold already applied)
zdna_regions = defaultdict(list)  # chrom -> [(start, end)]
with open(BEDGRAPH) as f:
    for line in f:
        if line.startswith("track") or line.startswith("browser"): continue
        p = line.split()
        if len(p) < 3: continue
        zdna_regions[p[0]].append((int(p[1]), int(p[2])))

print(f"  Z-DNA regions loaded: {sum(len(v) for v in zdna_regions.values())}")

# Build O.bimaculoides gene CDS coords
obi_gene_coords = {}  # prod_norm -> (chrom, tss, strand)
with gzip.open(GFF_DIR / "O_bimaculoides.gff.gz", "rt", errors="replace") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        if len(p) < 9 or p[2] != "CDS": continue
        prod_m = re.search(r'product=([^;]+)', p[8])
        if not prod_m: continue
        prod = normalize(prod_m.group(1))
        if prod not in obi_gene_coords:
            chrom = p[0]; start, end = int(p[3]), int(p[4]); strand = p[6]
            tss = start if strand == "+" else end
            obi_gene_coords[prod] = (chrom, tss, strand)

def has_zdna_in_promoter(chrom, tss, strand, up=2000, down=200):
    """Check bedgraph for Z-DNA in promoter window."""
    if strand == "+": prom_s = max(0, tss-up); prom_e = tss+down
    else: prom_s = max(0, tss-down); prom_e = tss+up
    for zs, ze in zdna_regions.get(chrom, []):
        if zs < prom_e and ze > prom_s:
            return True
    return False

# ── Step 5: Find candidates with Z-DNA in all 3 species ─────────────────────
print("\n=== Finding 5 genes with Z-DNA in all 3 species ===")

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

# Build O.vulgaris GFF index for fetching coords
print("  Building O.vulgaris GFF index...")
vul_idx = defaultdict(list)
seen = set()
with gzip.open(GFF_DIR / "O_vulgaris.gff.gz", "rt", errors="replace") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        if len(p) < 9 or p[2] != "CDS": continue
        prod_m = re.search(r'product=([^;]+)', p[8])
        if not prod_m: continue
        prod = normalize(prod_m.group(1))
        chrom = p[0]; start, end = int(p[3]), int(p[4]); strand = p[6]
        tss = start if strand == "+" else end
        key = (prod, chrom, tss)
        if key not in seen:
            seen.add(key)
            vul_idx[prod].append((chrom, tss, strand))

# Build A.fangsiao GFF index
print("  Building A.fangsiao GFF index...")
fang_idx = defaultdict(list)
seen = set()
with gzip.open(GFF_DIR / "A_fangsiao.gff.gz", "rt", errors="replace") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        if len(p) < 9 or p[2] != "CDS": continue
        prod_m = re.search(r'product=([^;]+)', p[8])
        if not prod_m: continue
        prod = normalize(prod_m.group(1))
        chrom = p[0]; start, end = int(p[3]), int(p[4]); strand = p[6]
        tss = start if strand == "+" else end
        key = (prod, chrom, tss)
        if key not in seen:
            seen.add(key)
            fang_idx[prod].append((chrom, tss, strand))

results = []
print(f"  Checking {len(all_candidates)} candidates...\n")

for i, cand in enumerate(all_candidates):
    # O. bimaculoides: check bedgraph
    obi_match = None
    for prod, coords in obi_gene_coords.items():
        if cand == prod or (len(cand) > 10 and (cand in prod or prod in cand)):
            if has_zdna_in_promoter(*coords):
                obi_match = (prod, coords)
                break
    if obi_match is None:
        if i < 5: print(f"  [{i}] {cand[:45]} — no Z-DNA in O.bimaculoides (bedgraph)")
        continue

    # O. vulgaris: find coords and fetch promoter
    vul_coords = None
    for vp, clist in vul_idx.items():
        if cand == vp or (len(cand) > 10 and (cand in vp or vp in cand)):
            vul_coords = clist[0]; break
    if vul_coords is None:
        if i < 5: print(f"  [{i}] {cand[:45]} — not in O.vulgaris GFF")
        continue

    seq_vul = fetch_promoter(*vul_coords)
    time.sleep(0.34)
    if seq_vul is None: continue
    zdna_vul = best_zdna(seq_vul)
    if zdna_vul is None:
        print(f"  [{i}] {cand[:45]} — no Z-DNA motif in O.vulgaris promoter"); continue

    # A. fangsiao: find coords and fetch promoter
    fang_coords = None
    for fp, clist in fang_idx.items():
        if cand == fp or (len(cand) > 10 and (cand in fp or fp in cand)):
            fang_coords = clist[0]; break
    if fang_coords is None:
        print(f"  [{i}] {cand[:45]} — not in A.fangsiao GFF"); continue

    seq_fang = fetch_promoter(*fang_coords)
    time.sleep(0.34)
    if seq_fang is None: continue
    zdna_fang = best_zdna(seq_fang)
    if zdna_fang is None:
        print(f"  [{i}] {cand[:45]} — no Z-DNA motif in A.fangsiao promoter"); continue

    # O. bimaculoides: fetch promoter sequence for the motif
    seq_obi = fetch_promoter(*obi_match[1])
    time.sleep(0.34)
    if seq_obi is None: continue
    zdna_obi = best_zdna(seq_obi)
    if zdna_obi is None: continue

    results.append({
        "gene": cand,
        "O_bimaculoides": zdna_obi,
        "O_vulgaris":     zdna_vul,
        "A_fangsiao":     zdna_fang,
    })
    print(f"  HIT [{len(results)}]: {cand}")
    if len(results) >= 5:
        print("  Got 5 genes, stopping!"); break

    if i % 5 == 0:
        print(f"  [{i}/{len(all_candidates)}] {len(results)} found so far...")

print(f"\n{'='*55}")
print(f"Genes with Z-DNA in all 3 species: {len(results)}")
for r in results:
    print(f"  {r['gene']}")
    for sp in ["O_bimaculoides","O_vulgaris","A_fangsiao"]:
        print(f"    {sp}: {r[sp]}")

rows = []
for r in results:
    for sp in ["O_bimaculoides","O_vulgaris","A_fangsiao"]:
        rows.append({"gene": r["gene"], "species": sp, "zdna_seq": r[sp]})

with open(BASE / "zdna_orthologs.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene","species","zdna_seq"], delimiter="\t")
    w.writeheader(); w.writerows(rows)
with open(BASE / "zdna_orthologs.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved zdna_orthologs.tsv ({len(rows)} rows)")
