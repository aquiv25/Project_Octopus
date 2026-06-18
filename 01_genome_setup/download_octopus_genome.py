"""
Download Octopus bimaculoides genome (GCF_001194135.2) from NCBI
"""

import os
import subprocess
import sys

# ----------------------------
# Config
# ----------------------------
ACCESSION = "GCF_001194135.2"
OUT_DIR = "octopus_genome"
GENOME_FILE = os.path.join(OUT_DIR, f"{ACCESSION}_ASM119413v2_genomic.fna.gz")
GENOME_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/194/135/"
    "GCF_001194135.2_ASM119413v2/"
    "GCF_001194135.2_ASM119413v2_genomic.fna.gz"
)

os.makedirs(OUT_DIR, exist_ok=True)


def run(cmd):
    print(f"$ {cmd}")
    ret = subprocess.run(cmd, shell=True)
    if ret.returncode != 0:
        sys.exit(f"Command failed: {cmd}")


# ----------------------------
# Option 1: wget (recommended)
# ----------------------------
def download_wget():
    run(f"wget -c -P {OUT_DIR} {GENOME_URL}")


# ----------------------------
# Option 2: NCBI datasets CLI (if installed)
# Gives checksums + structured metadata automatically
# Install: conda install -c conda-forge ncbi-datasets-cli
# ----------------------------
def download_datasets_cli():
    run(
        f"datasets download genome accession {ACCESSION} "
        f"--include genome --filename {OUT_DIR}/{ACCESSION}.zip"
    )
    run(f"unzip -o {OUT_DIR}/{ACCESSION}.zip -d {OUT_DIR}")


# ----------------------------
# Decompress
# ----------------------------
def decompress():
    gz = GENOME_FILE
    fna = gz.replace(".gz", "")
    if os.path.exists(gz) and not os.path.exists(fna):
        run(f"gunzip -k {gz}")
    return fna


if __name__ == "__main__":
    method = sys.argv[1] if len(sys.argv) > 1 else "wget"

    if method == "datasets":
        download_datasets_cli()
        fna_path = os.path.join(OUT_DIR, "ncbi_dataset/data", ACCESSION,
                                f"{ACCESSION}_ASM119413v2_genomic.fna")
    else:
        download_wget()
        fna_path = decompress()

    print(f"\nGenome ready at: {os.path.abspath(fna_path)}")
    print("\nUpdate your notebook:")
    print(f'  fasta_sequences = SeqIO.parse(open("{os.path.abspath(fna_path)}"), "fasta")')
