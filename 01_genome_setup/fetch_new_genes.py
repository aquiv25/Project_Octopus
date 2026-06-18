"""
Скачивает протеом Argonauta hians, извлекает 5 новых генов
из всех видов (включая уже скачанные протеомы),
сохраняет по одному FASTA на ген для MAFFT + IQ-TREE.

5 новых генов (из O. bimaculoides):
  LOC106870414 → Tropomyosin
  LOC106878893 → Mediator of RNA polymerase II transcription subunit 28
  LOC106872502 → SIK1 (serine/threonine-protein kinase SIK1)
  LOC106876456 → Homeobox protein abdominal-A
  LOC106882331 → Zinc finger protein 271
"""
import subprocess, gzip, shutil
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

BASE  = Path("/Users/dassagaripova/Documents/Claude/Projects/Article")
PROT  = BASE / "proteomes"
OUT   = BASE / "per_gene_new"
PROT.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)

# ── Уже скачанные протеомы (используем как есть) ────────────────────────
EXISTING_PROTEOMES = {
    "O_bimaculoides": PROT / "O_bimaculoides.faa",
    "O_sinensis":     PROT / "O_sinensis.faa",
    "A_fangsiao":     PROT / "A_fangsiao.faa",
    "O_rubescens":    PROT / "O_rubescens.faa",
    "O_vulgaris":     PROT / "O_vulgaris_full.faa",
    # E_cirrhosa: no full proteome available on NCBI FTP (genome-only assembly)
}

# ── Argonauta hians (GCA — не RefSeq, но содержит белки) ────────────────
NEW_SPECIES = {
    "A_hians": ("GCA_054771915.1", "GCA_054771915.1_Ahia01"),
}

# ── 5 новых генов ────────────────────────────────────────────────────────
TARGETS = {
    "TPM":    ["tropomyosin"],          # Tropomyosin
    "MED28":  ["mediator of rna polymerase ii transcription subunit 28",
               "mediator subunit 28", "med28"],
    "SIK1":   ["serine/threonine-protein kinase sik1", "salt inducible kinase 1",
               "sik1"],
    "ABDA":   ["homeobox protein abdominal-a", "abdominal-a", "abd-a"],
    "ZNF271": ["zinc finger protein 271", "znf271"],
}

# ════════════════════════════════════════════════════════════════════════
def download_proteome(sp_name, accession, known_folder):
    out_path = PROT / f"{sp_name}.faa"
    if out_path.exists():
        print(f"⏭  {sp_name} — already downloaded")
        return out_path

    # Handle both GCF and GCA
    prefix = accession[:3]   # GCF or GCA
    num = accession.replace(f"{prefix}_", "").split(".")[0]
    p1, p2, p3 = num[:3], num[3:6], num[6:9]
    ftp_base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{p1}/{p2}/{p3}/"

    faa_url = f"{ftp_base}{known_folder}/{known_folder}_protein.faa.gz"
    gz_path = PROT / f"{sp_name}.faa.gz"
    print(f"⬇  Downloading {sp_name} from {faa_url}")
    ret = subprocess.run(
        ['curl', '-s', '-L', faa_url, '-o', str(gz_path)],
        capture_output=True, timeout=180
    )
    if gz_path.exists() and gz_path.stat().st_size > 10_000:
        with gzip.open(gz_path, "rb") as fin, open(out_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        gz_path.unlink()
        n = int(subprocess.run(f'grep -c ">" "{out_path}"',
                               shell=True, capture_output=True, text=True).stdout.strip() or 0)
        print(f"✅  {sp_name}: {n:,} proteins")
        return out_path
    else:
        if gz_path.exists(): gz_path.unlink()
        print(f"❌  {sp_name}: download failed — {faa_url}")
        return None


def extract_genes(faa_path, sp_name):
    """Return dict gene -> SeqRecord for one species."""
    found = {}
    for record in SeqIO.parse(str(faa_path), "fasta"):
        # Check ID prefix (E. cirrhosa format: GENE|acc|species)
        id_prefix = record.id.split("|")[0].upper()
        if id_prefix in TARGETS and id_prefix not in found:
            rec = SeqRecord(record.seq, id=f"{id_prefix}_{sp_name}", description="")
            found[id_prefix] = rec
            continue
        # Keyword search in description
        text = (record.id + " " + record.description).lower()
        for gene, keywords in TARGETS.items():
            if gene not in found:
                if any(kw.lower() in text for kw in keywords):
                    rec = SeqRecord(record.seq, id=f"{gene}_{sp_name}", description="")
                    found[gene] = rec
    return found


# ════════════════════════════════════════════════════════════════════════
print("=" * 60)
print("STEP 1: Download A. hians proteome")
print("=" * 60)
for sp, (acc, folder) in NEW_SPECIES.items():
    download_proteome(sp, acc, folder)

print("\n" + "=" * 60)
print("STEP 2: Extract genes from all species")
print("=" * 60)

all_seqs = {g: {} for g in TARGETS}

# All proteomes to process
all_proteomes = {**EXISTING_PROTEOMES}
for sp in NEW_SPECIES:
    all_proteomes[sp] = PROT / f"{sp}.faa"

for sp_name, faa_path in all_proteomes.items():
    if not faa_path.exists():
        print(f"⚠️  {sp_name}: file not found ({faa_path})")
        continue
    found = extract_genes(faa_path, sp_name)
    for gene, rec in found.items():
        all_seqs[gene][sp_name] = rec
    print(f"{sp_name}: found {len(found)}/5 genes — {sorted(found.keys())}")

# ════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("STEP 3: Save per-gene FASTA files")
print("=" * 60)

SPECIES_ORDER = [
    "O_bimaculoides", "O_vulgaris", "O_sinensis",
    "A_fangsiao", "E_cirrhosa", "A_hians",
]

for gene, sp_dict in all_seqs.items():
    out_file = OUT / f"{gene}_species.fasta"
    records = []
    for sp in SPECIES_ORDER:
        if sp in sp_dict:
            records.append(sp_dict[sp])
        else:
            print(f"  ⚠️  {gene}: no sequence for {sp}")
    from Bio import SeqIO as _SeqIO
    _SeqIO.write(records, out_file, "fasta")
    print(f"📄  {gene}: {len(records)} species → {out_file.name}")

# Summary table
print("\n" + "=" * 60)
print("SUMMARY:")
print(f"{'Gene':<10}", end="")
for sp in SPECIES_ORDER:
    print(f"  {sp[:10]:<10}", end="")
print()
print("-" * (10 + 12 * len(SPECIES_ORDER)))
for gene in TARGETS:
    print(f"{gene:<10}", end="")
    for sp in SPECIES_ORDER:
        status = "✅" if sp in all_seqs[gene] else "❌"
        print(f"  {status:<10}", end="")
    print()

print(f"\n✅ Files saved to: {OUT}")
