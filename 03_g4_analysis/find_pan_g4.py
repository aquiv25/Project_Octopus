"""
Fast search: find 5 orthologous genes with strict G4 in promoters across all 5 species.
Uses threading for parallel NCBI fetches.
"""
import re, gzip, time, urllib.request, json, csv, sys
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.Seq import Seq

BASE = Path("/Users/dassagaripova/Documents/Claude/Projects/Article/g4_promoter")
G4_STRICT = re.compile(r'G{3,6}[ACGT]{1,7}G{3,6}[ACGT]{1,7}G{3,6}[ACGT]{1,7}G{3,6}', re.IGNORECASE)
ALL_SPECIES = ["O_bimaculoides","O_vulgaris","O_sinensis","A_fangsiao","A_hians"]

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

all_prods = [set(sp_idx[sp].keys()) for sp in ALL_SPECIES]
common = all_prods[0]
for s in all_prods[1:]: common &= s
named = sorted(p for p in common
               if "uncharacterized" not in p
               and "hypothetical" not in p
               and len(p) > 8)
print(f"Common named products: {len(named)}", flush=True)

def fetch_and_check(chrom, tss, strand):
    if strand == "+": s = max(1, tss-3000); e = tss+300
    else: s = max(1, tss-300); e = tss+3000
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={chrom}&seq_start={s}&seq_stop={e}"
           f"&rettype=fasta&retmode=text&email=researcher@example.com")
    try:
        with urllib.request.urlopen(url, timeout=12) as r:
            text = r.read().decode()
        if not text.startswith(">"): return None
        seq = "".join(text.strip().split("\n")[1:]).upper()
        if strand == "-": seq = str(Seq(seq).reverse_complement())
        m = G4_STRICT.search(seq) or G4_STRICT.search(str(Seq(seq).reverse_complement()))
        return m.group() if m else None
    except:
        return None

results = []
print("Scanning candidates...", flush=True)

for i, prod in enumerate(named):
    print(f"[{i}/{len(named)}] {prod[:50]}...", flush=True)
    sp_g4 = {}
    failed = False

    # fetch all 5 species in parallel
    with ThreadPoolExecutor(max_workers=5) as ex:
        futures = {}
        for sp in ALL_SPECIES:
            chrom, tss, strand = sp_idx[sp][prod][0]
            futures[ex.submit(fetch_and_check, chrom, tss, strand)] = sp
        for fut in as_completed(futures):
            sp = futures[fut]
            res = fut.result()
            if res is None:
                failed = True
            else:
                sp_g4[sp] = res

    if not failed and len(sp_g4) == 5:
        results.append({"product": prod, "g4": sp_g4})
        print(f"  >>> HIT: {prod}", flush=True)
        if len(results) >= 5:
            print("Got 5 genes, done!", flush=True)
            break
    else:
        print(f"  skip (failed={failed}, hits={len(sp_g4)}/5)", flush=True)

print(f"\nFinal: {len(results)} genes with strict G4 in all 5 species", flush=True)
for r in results:
    print(f"  {r['product']}", flush=True)

with open(BASE/"results/pan_g4_genes.json","w") as f:
    json.dump(results, f, indent=2)

rows = [{"gene": r["product"], "species": sp, "G4_seq": g4, "type": "strict"}
        for r in results for sp, g4 in r["g4"].items()]
with open(BASE/"results/pan_g4_strict.tsv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene","species","G4_seq","type"], delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"Saved {len(rows)} rows to pan_g4_strict.tsv", flush=True)
