"""
01_download_genomes.py
=======================
Downloads the SARS-CoV-2 and E. coli genomes from NCBI
using Biopython's Entrez API.

Saves:
    data/sars_cov2.fasta   — SARS-CoV-2 genome (NC_045512)
    data/ecoli.fasta       — E. coli K-12 genome (U00096)

Usage:
    python 01_download_genomes.py --email your.email@example.com

Author: Vedika Judyani
"""

import argparse
import os
from Bio import Entrez, SeqIO


def download_genome(accession, email, output_path):
    """
    Fetch a genome sequence from NCBI Nucleotide database.

    Parameters
    ----------
    accession   : str  — NCBI accession number (e.g. 'NC_045512')
    email       : str  — required by NCBI for API access
    output_path : str  — where to save the FASTA file
    """
    Entrez.email = email
    print(f"[INFO] Downloading {accession} from NCBI...")

    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="fasta",
        retmode="text"
    )
    sequence = handle.read()
    handle.close()

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write(sequence)

    print(f"[INFO] Saved to: {output_path}")
    return output_path


def verify_fasta(filepath):
    """
    Quick check — load the FASTA and print basic stats.
    """
    record = SeqIO.read(filepath, "fasta")
    print(f"[INFO] ID     : {record.id}")
    print(f"[INFO] Length : {len(record.seq):,} bp")
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download genomes from NCBI.")
    parser.add_argument("--email", required=True,
                        help="Your email address (required by NCBI Entrez)")
    args = parser.parse_args()

    # ── SARS-CoV-2 (RefSeq: NC_045512.2) ──
    sars_path = download_genome(
        accession="NC_045512",
        email=args.email,
        output_path="data/sars_cov2.fasta"
    )
    verify_fasta(sars_path)

    print()

    # ── E. coli K-12 MG1655 (RefSeq: U00096) ──
    ecoli_path = download_genome(
        accession="U00096",
        email=args.email,
        output_path="data/ecoli.fasta"
    )
    verify_fasta(ecoli_path)

    print("\n[DONE] Both genomes downloaded successfully.")
    print("       Next step: run  python 02_orf_prediction.py")
