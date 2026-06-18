"""
Find G4 in promoters of 5 target genes across all 5 species.
Uses http (not https) to avoid SSL issues, with retry logic.
Falls back to relaxed G4 if strict not found.
"""
import re, gzip, time, urllib.request, json, csv
from pathlib import Path
from collections import defaultdict
from Bio.Seq import Seq

BASE = Path("/Users/dassagaripova/Documents/Claude/Projects/Article/g4_promoter")
G4_STRICT  = re.compile(r'G{3,6}[ACGT]{1,7}G{3,6}[ACGT]{1,7}G{3,6}[ACGT]{1,7}G{3,6}', re.IGNORECASE)
G4_RELAXED = re.compile(r'G{2,6}[ACGT]{1,12}G{2,6}[ACGT]{1,12}G{2,6}[ACGT]{1,12}G{2,6}', re.IGNORECASE)
ALL_SPECIES = ["O_bimaculoides","O_vulgaris","O_sinensis","A_fangsiao","A_hians"]

TARGETS = {
    "TPM":    ["tropomyosin"],
    "MED28":  ["mediator of rna polymerase ii transcription subuni", "mediator complex subunit 28", "med28"],
    "SIK1":   ["serine/threonine-protein kinase sik1", "salt inducible kinase", "sik1"],
    "ABDA":   ["homeobox protein abdominal-a", "homeobox abdominal", "abdominal-a", "abd-a",
               "homeobox protein abd", "abdominal a"],
    "ZNF271": ["zinc finger protein 271"],
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

def find_gene(kws, idx):
    for prod, coords in idx.items():
        for kw in kws:
            if kw in prod:
                return prod, coords[0]
    return None, None

def fetch_promoter(chrom, tss, strand, up=3000, down=300):
    if strand == "+": s = max(1, tss-up); e = tss+down
    else: s = max(1, tss-down); e = tss+up
    # Use http to avoid SSL issues
    url = (f"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={chrom}&seq_start={s}&seq_stop={e}"
           f"&rettype=fasta&retmode=text&email=researcher@example.com")
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                text = r.read().decode()
            if not text.startswith(">"): return None
            seq = "".join(text.strip().split("\n")[1:]).upper()
            if strand == "-": seq = str(Seq(seq).reverse_complement())
            return seq
        except Exception as e:
            print(f"    attempt {attempt+1} failed: {e}", flush=True)
            time.sleep(2)
    return None

def find_g4(seq):
    """Returns (g4_seq, type) or (None, None)"""
    rc = str(Seq(seq).reverse_complement())
    for pattern, ptype in [(G4_STRICT, "strict"), (G4_RELAXED, "relaxed")]:
        hits = list(pattern.finditer(seq)) + list(pattern.finditer(rc))
        if hits:
            best = max(hits, key=lambda m: len(m.group()))
            return best.group(), ptype
    return None, None

# ── Main ─────────────────────────────────────────────────────────────────────
print("\nFetching promoters and searching G4...", flush=True)
rows = []

for gene, kws in TARGETS.items():
    print(f"\n=== {gene} ===", flush=True)
    for sp in ALL_SPECIES:
        prod_name, coords = find_gene(kws, sp_idx[sp])
        if coords is None:
            print(f"  {sp}: NOT FOUND", flush=True)
            continue
        chrom, tss, strand = coords
        print(f"  {sp}: {prod_name[:55]} | {chrom} {strand}", flush=True)
        seq = fetch_promoter(chrom, tss, strand)
        time.sleep(0.34)
        if seq is None:
            print(f"    -> fetch FAILED", flush=True)
            continue
        g4, gtype = find_g4(seq)
        if g4:
            print(f"    -> {gtype} G4: {g4}", flush=True)
            rows.append({"gene": gene, "species": sp, "G4_seq": g4,
                         "type": gtype, "product": prod_name})
        else:
            print(f"    -> no G4 found", flush=True)

print(f"\n{'='*55}", flush=True)
print(f"Total: {len(rows)} hits", flush=True)
for r in rows:
    print(f"  {r['gene']}/{r['species']} [{r['type']}]: {r['G4_seq']}", flush=True)

with open(BASE/"results/target_g4_v2.tsv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene","species","G4_seq","type","product"], delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"Saved results/target_g4_v2.tsv", flush=True)
