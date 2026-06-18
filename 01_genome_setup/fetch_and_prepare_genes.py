"""
Скачивает протеомы 6 видов с NCBI, извлекает 8 эпигенетических генов,
объединяет с уже имеющимися O. vulgaris + E. cirrhosa,
сохраняет по одному FASTA на ген — готово для MEGA.
"""
import os, subprocess, gzip, shutil
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ── Пути ────────────────────────────────────────────────────────────────
BASE  = Path("/Users/dassagaripova/Documents/Claude/Projects/Article")
PROT  = BASE / "proteomes"
OUT   = BASE / "per_gene_mega"
PROT.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)

# ── Уже имеющиеся файлы ─────────────────────────────────────────────────
EXISTING = {
    "O_vulgaris":  Path("/var/folders/1l/9_k40zy95y7_ymxby1kwr6v40000gn/T/octopus_vulgaris_8genes.faa"),
    "E_cirrhosa":  Path("/var/folders/1l/9_k40zy95y7_ymxby1kwr6v40000gn/T/epigenetic_8genes.fasta"),
}

# ── 6 видов для скачивания ───────────────────────────────────────────────
# accession → точное имя папки на FTP (чтобы не делать листинг)
SPECIES_TO_DOWNLOAD = {
    "O_bimaculoides": ("GCF_001194135.2", "GCF_001194135.2_ASM119413v2"),
    "O_americanus":   ("GCF_030220305.1", None),
    "O_sinensis":     ("GCF_006345805.1", "GCF_006345805.1_ASM634580v1"),
    "A_fangsiao":     ("GCF_014858855.1", "GCF_014858855.1_ASM1485885v1"),
    "O_rubescens":    ("GCF_025177155.1", "GCF_025177155.1_ASM2517715v1"),  # исправлен аккессион
    "O_mimus":        ("GCF_964030095.1", None),
}

# ── Гены и ключевые слова для поиска ────────────────────────────────────
TARGETS = {
    "DNMT1": ["cytosine-5)-methyltransferase 1", "dna methyltransferase 1"],
    "TET3":  ["methylcytosine dioxygenase tet3", "tet3", "methylcytosine dioxygenase TET"],
    "EZH2":  ["enhancer of zeste homolog 2", "histone-lysine n-methyltransferase ezh2",
               "histone-lysine N-methyltransferase EZH2"],
    "KMT2A": ["lysine n-methyltransferase 2a", "histone-lysine n-methyltransferase 2a",
               "myeloid/lymphoid or mixed-lineage leukemia"],
    "CHD1":  ["chromodomain-helicase-dna-binding protein 1", "chd1"],
    "BRD3":  ["bromodomain-containing protein 3"],
    "ADAR1": ["double-stranded rna-specific adenosine deaminase",
               "rna-specific adenosine deaminase", "adenosine deaminase adar"],
    "UHRF1": ["ubiquitin-protein ligase uhrf1", "e3 ubiquitin-protein ligase uhrf1",
               "uhrf1"],
}

# ════════════════════════════════════════════════════════════════════════
# 1. Скачать протеомы
# ════════════════════════════════════════════════════════════════════════
def download_proteome(sp_name, accession, known_folder=None):
    out_path = PROT / f"{sp_name}.faa"
    if out_path.exists():
        print(f"⏭  {sp_name} — уже скачан ({out_path})")
        return out_path

    num = accession.replace("GCF_", "").split(".")[0]
    p1, p2, p3 = num[:3], num[3:6], num[6:9]
    ftp_base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{p1}/{p2}/{p3}/"

    import re
    folder = known_folder or ""

    if not folder:
        # Найти папку через HTML листинг FTP (с увеличенным таймаутом)
        try:
            result = subprocess.run(
                ['curl', '-s', '--max-time', '60', ftp_base],
                capture_output=True, text=True, timeout=65
            )
            matches = re.findall(rf'{re.escape(accession)}[^"</ \n]*', result.stdout)
            for m in matches:
                m = m.rstrip("/")
                if "_" in m[len(accession):]:
                    folder = m
                    break
        except Exception as e:
            print(f"  ⚠️  FTP листинг не удался: {e}")

    if not folder:
        print(f"❌  {sp_name}: папка не найдена на FTP")
        return None

    faa_url = f"{ftp_base}{folder}/{folder}_protein.faa.gz"
    gz_path = PROT / f"{sp_name}.faa.gz"
    ret = subprocess.run(
        ['curl', '-s', '-L', faa_url, '-o', str(gz_path)],
        capture_output=True, timeout=120
    )
    if gz_path.exists() and gz_path.stat().st_size > 10_000:
        with gzip.open(gz_path, "rb") as fin, open(out_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        gz_path.unlink()
        n = int(subprocess.run(f'grep -c ">" "{out_path}"',
                               shell=True, capture_output=True, text=True).stdout.strip() or 0)
        print(f"✅  {sp_name}: {n:,} белков скачано")
        return out_path
    else:
        if gz_path.exists(): gz_path.unlink()
        print(f"❌  {sp_name}: ошибка скачивания — {faa_url}")
        return None

print("=" * 60)
print("ШАГ 1: Скачиваем протеомы")
print("=" * 60)
for sp, (acc, known_folder) in SPECIES_TO_DOWNLOAD.items():
    download_proteome(sp, acc, known_folder)

# ════════════════════════════════════════════════════════════════════════
# 2. Извлечь гены из скачанных протеомов
# ════════════════════════════════════════════════════════════════════════
def extract_genes(faa_path, sp_name):
    """Вернуть dict gene -> SeqRecord для одного вида."""
    found = {}
    for record in SeqIO.parse(faa_path, "fasta"):
        # Сначала проверить: ID начинается с имени гена (формат E. cirrhosa)
        # Пример: EZH2|ENSHJAP00000004189.1|Eledone_cirrhosa
        id_prefix = record.id.split("|")[0].upper()
        if id_prefix in TARGETS and id_prefix not in found:
            rec = SeqRecord(record.seq, id=f"{id_prefix}_{sp_name}", description="")
            found[id_prefix] = rec
            continue

        # Иначе — поиск по ключевым словам в описании
        text = (record.id + " " + record.description).lower()
        for gene, keywords in TARGETS.items():
            if gene not in found:
                if any(kw.lower() in text for kw in keywords):
                    rec = SeqRecord(record.seq, id=f"{gene}_{sp_name}", description="")
                    found[gene] = rec
    return found

print("\n" + "=" * 60)
print("ШАГ 2: Извлекаем гены")
print("=" * 60)

# gene -> {sp_name: SeqRecord}
all_seqs = {g: {} for g in TARGETS}

# --- Обработать уже имеющиеся файлы ---
for sp_name, faa_path in EXISTING.items():
    if not faa_path.exists():
        print(f"⚠️  Файл не найден: {faa_path}")
        continue
    found = extract_genes(faa_path, sp_name)
    for gene, rec in found.items():
        all_seqs[gene][sp_name] = rec
    print(f"{sp_name}: найдено {len(found)}/8 генов — {sorted(found.keys())}")

# --- Обработать скачанные протеомы ---
for sp_name in SPECIES_TO_DOWNLOAD:
    faa_path = PROT / f"{sp_name}.faa"
    if not faa_path.exists():
        print(f"⚠️  {sp_name}: файл не скачан, пропускаю")
        continue
    found = extract_genes(faa_path, sp_name)
    for gene, rec in found.items():
        all_seqs[gene][sp_name] = rec
    print(f"{sp_name}: найдено {len(found)}/8 генов — {sorted(found.keys())}")

# ════════════════════════════════════════════════════════════════════════
# 3. Сохранить по одному FASTA на ген
# ════════════════════════════════════════════════════════════════════════
SPECIES_ORDER = [
    "O_bimaculoides", "O_americanus", "O_vulgaris", "O_sinensis",
    "A_fangsiao", "O_rubescens", "O_mimus", "E_cirrhosa"
]

print("\n" + "=" * 60)
print("ШАГ 3: Сохраняем файлы для MEGA")
print("=" * 60)

for gene, sp_dict in all_seqs.items():
    out_file = OUT / f"{gene}_8species.fasta"
    records = []
    for sp in SPECIES_ORDER:
        if sp in sp_dict:
            records.append(sp_dict[sp])
        else:
            print(f"  ⚠️  {gene}: нет последовательности для {sp}")
    SeqIO.write(records, out_file, "fasta")
    print(f"📄  {gene}: {len(records)}/8 видов → {out_file.name}")

# ── Сводная таблица ─────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("ИТОГОВАЯ ТАБЛИЦА:")
print(f"{'Ген':<10}", end="")
for sp in SPECIES_ORDER:
    print(f"  {sp[:8]:<8}", end="")
print()
print("-" * (10 + 10 * len(SPECIES_ORDER)))
for gene in TARGETS:
    print(f"{gene:<10}", end="")
    for sp in SPECIES_ORDER:
        status = "✅" if sp in all_seqs[gene] else "❌"
        print(f"  {status:<8}", end="")
    print()

print(f"\n✅ Файлы сохранены в: {OUT}")
print("📌 Следующий шаг: открыть каждый файл в MEGA → Align → Align by MUSCLE → строить дерево")
