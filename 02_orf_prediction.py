"""
02_orf_prediction.py
=====================
Identifies Open Reading Frames (ORFs) in a genome sequence.

Pipeline:
    1. First validates the ORF finder on E. coli (well-annotated ground truth)
    2. Then applies it to the SARS-CoV-2 genome
    3. Saves results as CSV and protein FASTA files

An ORF is a stretch of DNA starting with ATG (start codon)
and ending with TAA, TAG, or TGA (stop codons), in any of
the 6 reading frames (3 forward, 3 reverse complement).

Usage:
    python 02_orf_prediction.py --genome data/sars_cov2.fasta \
                                 --organism sars_cov2 \
                                 --min_length 100

Author: Vedika Judyani
"""

import argparse
import os
import csv
from Bio import SeqIO
from Bio.Seq import Seq


# ─────────────────────────────────────────────
# CODON TABLES
# ─────────────────────────────────────────────

START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}


# ─────────────────────────────────────────────
# 1. ORF DETECTION
# ─────────────────────────────────────────────

def find_orfs(seq_record, min_length=100):
    """
    Identify all ORFs in all 6 reading frames of a sequence.

    Reading frames searched:
        Forward : +1, +2, +3
        Reverse : -1, -2, -3 (reverse complement)

    Parameters
    ----------
    seq_record : Bio.SeqRecord
    min_length : int — minimum ORF length in nucleotides (default 100)

    Returns
    -------
    orfs : list of dicts with keys:
        orf_id, start, end, strand, frame, length_nt,
        length_aa, dna_seq, protein_seq
    """
    orfs = []
    seq = seq_record.seq.upper()
    seq_len = len(seq)

    # Search both strands
    strands = [
        (seq, "+"),
        (seq.reverse_complement(), "-")
    ]

    for strand_seq, strand in strands:
        for frame in range(3):
            i = frame
            while i < len(strand_seq) - 2:
                codon = str(strand_seq[i:i+3])

                if codon == START_CODON:
                    # Found start — scan for stop codon
                    start_pos = i
                    for j in range(i + 3, len(strand_seq) - 2, 3):
                        stop_codon = str(strand_seq[j:j+3])
                        if stop_codon in STOP_CODONS:
                            orf_dna = strand_seq[start_pos:j+3]
                            length_nt = len(orf_dna)

                            if length_nt >= min_length:
                                protein = Seq(str(orf_dna)).translate(to_stop=True)

                                # Compute genomic coordinates
                                if strand == "+":
                                    genomic_start = start_pos
                                    genomic_end   = j + 3
                                else:
                                    genomic_start = seq_len - (j + 3)
                                    genomic_end   = seq_len - start_pos

                                orfs.append({
                                    "orf_id":      f"ORF_{len(orfs)+1:05d}",
                                    "start":        genomic_start,
                                    "end":          genomic_end,
                                    "strand":       strand,
                                    "frame":        frame + 1,
                                    "length_nt":    length_nt,
                                    "length_aa":    len(protein),
                                    "dna_seq":      str(orf_dna),
                                    "protein_seq":  str(protein)
                                })
                            i = j + 3  # skip past this ORF
                            break
                    else:
                        i += 3
                else:
                    i += 3

    return orfs


# ─────────────────────────────────────────────
# 2. FILTER & EXPORT
# ─────────────────────────────────────────────

def filter_orfs(orfs, min_aa=33):
    """
    Remove ORFs that are too short to encode functional proteins.
    Standard threshold: ≥100 nt (≥33 amino acids).

    Parameters
    ----------
    orfs   : list of dicts
    min_aa : int — minimum protein length in amino acids

    Returns
    -------
    filtered : list of dicts
    """
    filtered = [o for o in orfs if o["length_aa"] >= min_aa]
    print(f"[INFO] ORFs before filter : {len(orfs):,}")
    print(f"[INFO] ORFs after filter  : {len(filtered):,}  (≥{min_aa} aa)")
    return filtered


def export_orf_table(orfs, output_csv):
    """
    Save ORF metadata to a CSV file (excludes raw sequences).

    Parameters
    ----------
    orfs       : list of dicts
    output_csv : str — output file path
    """
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    fields = ["orf_id", "start", "end", "strand", "frame",
              "length_nt", "length_aa"]

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for orf in orfs:
            writer.writerow({k: orf[k] for k in fields})

    print(f"[INFO] ORF table saved  : {output_csv}")


def export_protein_fasta(orfs, output_fasta):
    """
    Export all ORF protein sequences to a FASTA file.

    Parameters
    ----------
    orfs         : list of dicts
    output_fasta : str — output file path
    """
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    with open(output_fasta, "w") as f:
        for orf in orfs:
            f.write(f">{orf['orf_id']} "
                    f"start={orf['start']} end={orf['end']} "
                    f"strand={orf['strand']} frame={orf['frame']} "
                    f"len={orf['length_aa']}aa\n")
            # Write sequence in 60-character lines (standard FASTA)
            seq = orf["protein_seq"]
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

    print(f"[INFO] Protein FASTA saved: {output_fasta}")


# ─────────────────────────────────────────────
# 3. VALIDATION ON E. COLI
# ─────────────────────────────────────────────

def validate_on_ecoli(ecoli_fasta, min_length=100):
    """
    Run ORF prediction on E. coli as a sanity check.
    E. coli K-12 has ~4,300 annotated protein-coding genes.
    We expect to find a comparable number of ORF candidates.

    Parameters
    ----------
    ecoli_fasta : str  — path to E. coli FASTA
    min_length  : int

    Returns
    -------
    orfs : list of dicts
    """
    print("\n" + "="*50)
    print("  VALIDATION: E. coli ORF Prediction")
    print("="*50)

    record = SeqIO.read(ecoli_fasta, "fasta")
    print(f"[INFO] E. coli genome length: {len(record.seq):,} bp")

    orfs = find_orfs(record, min_length=min_length)
    orfs = filter_orfs(orfs)

    print(f"[INFO] E. coli has ~4,300 annotated genes.")
    print(f"[INFO] Our prediction found {len(orfs):,} ORF candidates.")
    print(f"[INFO] False positives expected — ORF finder is intentionally")
    print(f"       permissive; downstream filtering reduces false positives.")

    export_orf_table(orfs, "results/ecoli_orfs.csv")
    export_protein_fasta(orfs, "results/ecoli_proteome.fasta")

    return orfs


# ─────────────────────────────────────────────
# 4. MAIN
# ─────────────────────────────────────────────

def main(genome_fasta, organism, min_length):
    print("\n" + "="*50)
    print(f"  ORF PREDICTION: {organism.upper()}")
    print("="*50)

    record = SeqIO.read(genome_fasta, "fasta")
    print(f"[INFO] Genome ID     : {record.id}")
    print(f"[INFO] Genome length : {len(record.seq):,} bp")

    orfs = find_orfs(record, min_length=min_length)
    orfs = filter_orfs(orfs)

    csv_out   = f"results/{organism}_orfs.csv"
    fasta_out = f"results/{organism}_proteome.fasta"

    export_orf_table(orfs, csv_out)
    export_protein_fasta(orfs, fasta_out)

    print(f"\n[DONE] {len(orfs):,} ORFs identified for {organism}")
    print(f"       Next step: run  python 03_hydrophobicity.py")
    return orfs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ORF prediction pipeline.")
    parser.add_argument("--genome",     required=True,
                        help="Path to genome FASTA file")
    parser.add_argument("--organism",   default="sars_cov2",
                        help="Organism name (used in output filenames)")
    parser.add_argument("--min_length", type=int, default=100,
                        help="Minimum ORF length in nucleotides (default: 100)")
    parser.add_argument("--validate",   action="store_true",
                        help="Also run validation on E. coli first")
    args = parser.parse_args()

    os.makedirs("results", exist_ok=True)

    if args.validate:
        validate_on_ecoli("data/ecoli.fasta", min_length=args.min_length)
        print()

    main(args.genome, args.organism, args.min_length)
