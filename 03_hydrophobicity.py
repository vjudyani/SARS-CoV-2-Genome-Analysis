"""
03_hydrophobicity.py
=====================
Computes hydrophobicity scores for all predicted ORF proteins
using the Kyte-Doolittle hydrophobicity scale.

High hydrophobicity → membrane-associated or secreted protein
Low hydrophobicity  → cytoplasmic / nuclear protein

For SARS-CoV-2, this helps identify:
  - Spike protein (membrane fusion, vaccine target)
  - Envelope protein (membrane protein)
  - Membrane protein (structural)

Outputs:
    results/sars_cov2_hydrophobicity.csv
    results/hydrophobicity_distribution.png
    results/top_hydrophobic_proteins.png

Usage:
    python 03_hydrophobicity.py \
        --fasta results/sars_cov2_proteome.fasta \
        --organism sars_cov2

Author: Vedika Judyani
"""

import argparse
import os
import csv
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


# ─────────────────────────────────────────────
# KYTE-DOOLITTLE HYDROPHOBICITY SCALE
# ─────────────────────────────────────────────
# Values reflect the tendency of each amino acid
# to be buried in a hydrophobic protein core.
# Source: Kyte & Doolittle (1982) J Mol Biol 157:105-132

KYTE_DOOLITTLE = {
    'A':  1.8,  # Alanine
    'R': -4.5,  # Arginine
    'N': -3.5,  # Asparagine
    'D': -3.5,  # Aspartic acid
    'C':  2.5,  # Cysteine
    'Q': -3.5,  # Glutamine
    'E': -3.5,  # Glutamic acid
    'G': -0.4,  # Glycine
    'H': -3.2,  # Histidine
    'I':  4.5,  # Isoleucine
    'L':  3.8,  # Leucine
    'K': -3.9,  # Lysine
    'M':  1.9,  # Methionine
    'F':  2.8,  # Phenylalanine
    'P': -1.6,  # Proline
    'S': -0.8,  # Serine
    'T': -0.7,  # Threonine
    'W': -0.9,  # Tryptophan
    'Y': -1.3,  # Tyrosine
    'V':  4.2,  # Valine
}


# ─────────────────────────────────────────────
# 1. PARSE PROTEIN FASTA
# ─────────────────────────────────────────────

def parse_protein_fasta(fasta_path):
    """
    Parse a protein FASTA file into a list of dicts.

    Parameters
    ----------
    fasta_path : str

    Returns
    -------
    proteins : list of dicts {orf_id, header, sequence}
    """
    proteins = []
    current_id = None
    current_header = None
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    proteins.append({
                        "orf_id":   current_id,
                        "header":   current_header,
                        "sequence": "".join(current_seq)
                    })
                parts = line[1:].split(" ", 1)
                current_id     = parts[0]
                current_header = line[1:]
                current_seq    = []
            else:
                current_seq.append(line)

    if current_id:
        proteins.append({
            "orf_id":   current_id,
            "header":   current_header,
            "sequence": "".join(current_seq)
        })

    print(f"[INFO] Loaded {len(proteins):,} protein sequences from {fasta_path}")
    return proteins


# ─────────────────────────────────────────────
# 2. COMPUTE HYDROPHOBICITY
# ─────────────────────────────────────────────

def compute_hydrophobicity(sequence, scale=KYTE_DOOLITTLE):
    """
    Compute mean hydrophobicity of a protein sequence
    using the Kyte-Doolittle scale.

    Unknown/non-standard amino acids are skipped.

    Parameters
    ----------
    sequence : str — single-letter amino acid sequence
    scale    : dict — hydrophobicity values per amino acid

    Returns
    -------
    mean_hydrophobicity : float
    """
    scores = [scale[aa] for aa in sequence.upper() if aa in scale]
    if not scores:
        return 0.0
    return round(sum(scores) / len(scores), 4)


def compute_window_hydrophobicity(sequence, window=9, scale=KYTE_DOOLITTLE):
    """
    Compute hydrophobicity along a sliding window.
    Used to identify transmembrane helices (typically 19-21 aa).

    Parameters
    ----------
    sequence : str
    window   : int — window size in amino acids
    scale    : dict

    Returns
    -------
    positions : list of int
    scores    : list of float
    """
    positions = []
    scores = []
    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i+window]
        score = compute_hydrophobicity(window_seq, scale)
        positions.append(i + window // 2)
        scores.append(score)
    return positions, scores


def analyze_proteins(proteins):
    """
    Compute hydrophobicity stats for all proteins.

    Parameters
    ----------
    proteins : list of dicts

    Returns
    -------
    df : pd.DataFrame with hydrophobicity metrics
    """
    results = []
    for p in proteins:
        seq = p["sequence"]
        mean_hydro = compute_hydrophobicity(seq)
        length_aa  = len(seq)

        # Classify by hydrophobicity
        if mean_hydro > 1.0:
            category = "Highly Hydrophobic (membrane-associated)"
        elif mean_hydro > 0.0:
            category = "Moderately Hydrophobic"
        elif mean_hydro > -1.0:
            category = "Hydrophilic"
        else:
            category = "Strongly Hydrophilic (cytoplasmic)"

        results.append({
            "orf_id":            p["orf_id"],
            "length_aa":         length_aa,
            "mean_hydrophobicity": mean_hydro,
            "category":          category
        })

    df = pd.DataFrame(results)
    return df


# ─────────────────────────────────────────────
# 3. VISUALIZATION
# ─────────────────────────────────────────────

def plot_hydrophobicity_distribution(df, organism, output_dir):
    """
    Plot distribution of mean hydrophobicity scores across all ORFs.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # ── Left: Histogram ──
    axes[0].hist(df["mean_hydrophobicity"], bins=40,
                 color='steelblue', edgecolor='white', alpha=0.85)
    axes[0].axvline(x=0, color='red', linestyle='--', linewidth=1.2,
                    label='Hydrophobic/Hydrophilic boundary')
    axes[0].axvline(x=df["mean_hydrophobicity"].mean(),
                    color='orange', linestyle='-', linewidth=1.5,
                    label=f'Mean = {df["mean_hydrophobicity"].mean():.3f}')
    axes[0].set_xlabel("Mean Kyte-Doolittle Hydrophobicity Score", fontsize=11)
    axes[0].set_ylabel("Number of ORFs", fontsize=11)
    axes[0].set_title(f"Hydrophobicity Distribution\n{organism.replace('_', ' ').title()}",
                      fontsize=13, fontweight='bold')
    axes[0].legend(fontsize=9)

    # ── Right: Category breakdown (pie chart) ──
    cat_counts = df["category"].value_counts()
    colors = ['#d73027', '#fc8d59', '#91bfdb', '#4575b4']
    wedges, texts, autotexts = axes[1].pie(
        cat_counts.values,
        labels=None,
        colors=colors[:len(cat_counts)],
        autopct='%1.1f%%',
        startangle=140,
        pctdistance=0.75
    )
    axes[1].set_title("ORF Category Breakdown\nby Hydrophobicity",
                      fontsize=13, fontweight='bold')
    legend_patches = [
        mpatches.Patch(color=colors[i], label=cat_counts.index[i])
        for i in range(len(cat_counts))
    ]
    axes[1].legend(handles=legend_patches, loc='lower center',
                   bbox_to_anchor=(0.5, -0.25), fontsize=8)

    plt.tight_layout()
    out_path = os.path.join(output_dir, f"{organism}_hydrophobicity_distribution.png")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Distribution plot saved: {out_path}")


def plot_top_hydrophobic(df, organism, output_dir, top_n=15):
    """
    Bar chart of top N most hydrophobic ORFs — these are
    candidates for membrane proteins and drug targets.
    """
    top = df.nlargest(top_n, "mean_hydrophobicity").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.barh(top["orf_id"], top["mean_hydrophobicity"],
                   color='#d73027', edgecolor='white', alpha=0.85)
    ax.axvline(x=0, color='black', linewidth=0.8)
    ax.set_xlabel("Mean Kyte-Doolittle Hydrophobicity Score", fontsize=11)
    ax.set_title(f"Top {top_n} Most Hydrophobic ORFs\n"
                 f"{organism.replace('_', ' ').title()} — Potential Membrane/Drug Target Candidates",
                 fontsize=12, fontweight='bold')

    # Annotate bars with length
    for i, row in top.iterrows():
        ax.text(row["mean_hydrophobicity"] + 0.02, i,
                f"{row['length_aa']} aa", va='center', fontsize=8)

    ax.invert_yaxis()
    plt.tight_layout()
    out_path = os.path.join(output_dir, f"{organism}_top_hydrophobic_orfs.png")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Top hydrophobic plot saved: {out_path}")


def plot_length_vs_hydrophobicity(df, organism, output_dir):
    """
    Scatter plot: protein length vs hydrophobicity score.
    Helps identify outliers (long + hydrophobic = likely membrane protein).
    """
    fig, ax = plt.subplots(figsize=(9, 6))

    scatter = ax.scatter(
        df["length_aa"],
        df["mean_hydrophobicity"],
        c=df["mean_hydrophobicity"],
        cmap='RdYlBu_r',
        alpha=0.6,
        s=25,
        edgecolors='none'
    )
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
    ax.set_xlabel("Protein Length (amino acids)", fontsize=11)
    ax.set_ylabel("Mean Hydrophobicity Score", fontsize=11)
    ax.set_title(f"Protein Length vs. Hydrophobicity\n"
                 f"{organism.replace('_', ' ').title()}",
                 fontsize=13, fontweight='bold')
    plt.colorbar(scatter, ax=ax, label='Hydrophobicity Score')
    plt.tight_layout()
    out_path = os.path.join(output_dir, f"{organism}_length_vs_hydrophobicity.png")
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[INFO] Scatter plot saved: {out_path}")


# ─────────────────────────────────────────────
# 4. MAIN
# ─────────────────────────────────────────────

def main(fasta_path, organism, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    print("\n" + "="*50)
    print(f"  HYDROPHOBICITY ANALYSIS: {organism.upper()}")
    print("="*50)

    # Parse proteins
    proteins = parse_protein_fasta(fasta_path)

    # Compute hydrophobicity
    df = analyze_proteins(proteins)

    # Print summary stats
    print(f"\n[INFO] Hydrophobicity Summary:")
    print(f"       Mean  : {df['mean_hydrophobicity'].mean():.3f}")
    print(f"       Std   : {df['mean_hydrophobicity'].std():.3f}")
    print(f"       Min   : {df['mean_hydrophobicity'].min():.3f}")
    print(f"       Max   : {df['mean_hydrophobicity'].max():.3f}")
    print(f"\n[INFO] Category breakdown:")
    print(df["category"].value_counts().to_string())

    # Save CSV
    csv_out = os.path.join(output_dir, f"{organism}_hydrophobicity.csv")
    df.to_csv(csv_out, index=False)
    print(f"\n[INFO] Results saved: {csv_out}")

    # Plots
    plot_hydrophobicity_distribution(df, organism, output_dir)
    plot_top_hydrophobic(df, organism, output_dir)
    plot_length_vs_hydrophobicity(df, organism, output_dir)

    print(f"\n[DONE] Hydrophobicity analysis complete.")
    print(f"       Next step: run  python 04_comparative_analysis.py")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Protein hydrophobicity analysis.")
    parser.add_argument("--fasta",    required=True,
                        help="Path to protein FASTA file (from ORF prediction)")
    parser.add_argument("--organism", default="sars_cov2",
                        help="Organism name for output labels")
    parser.add_argument("--output",   default="results/",
                        help="Output directory for plots and CSV")
    args = parser.parse_args()

    main(args.fasta, args.organism, args.output)
