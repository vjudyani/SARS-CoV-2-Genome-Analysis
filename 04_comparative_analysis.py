"""
04_comparative_analysis.py
===========================
Comparative analysis of SARS-CoV-2 against other viruses
using NCBI Virus metadata.

Generates:
    - Genome size distribution histogram by host category
    - Summary statistics table
    - GC content comparison across virus families

Input:
    data/viral_metadata.csv  — downloaded from NCBI Virus portal
    (Instructions in README for how to download this file)

Usage:
    python 04_comparative_analysis.py \
        --metadata data/viral_metadata.csv

    # Or use built-in synthetic data for testing:
    python 04_comparative_analysis.py --use_synthetic

Author: Vedika Judyani
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")


# ─────────────────────────────────────────────
# 1. GENERATE SYNTHETIC METADATA (for testing)
# ─────────────────────────────────────────────

def generate_synthetic_metadata(n_per_host=20, random_state=42):
    """
    Generate a synthetic viral metadata table that mimics the
    NCBI Virus download format.

    Host categories: Human, Bat, Avian, Rodent, Other mammal
    Families: Coronaviridae, Flaviviridae, Orthomyxoviridae,
              Retroviridae, Herpesviridae

    Returns
    -------
    df : pd.DataFrame
    """
    np.random.seed(random_state)

    host_categories = {
        "Human":        {"size_mean": 15000, "size_std": 12000,
                         "gc_mean": 0.48, "gc_std": 0.06},
        "Bat":          {"size_mean": 22000, "size_std": 8000,
                         "gc_mean": 0.40, "gc_std": 0.05},
        "Avian":        {"size_mean": 12000, "size_std": 5000,
                         "gc_mean": 0.50, "gc_std": 0.04},
        "Rodent":       {"size_mean": 10000, "size_std": 4000,
                         "gc_mean": 0.45, "gc_std": 0.05},
        "Other mammal": {"size_mean": 18000, "size_std": 9000,
                         "gc_mean": 0.44, "gc_std": 0.06},
    }

    families = ["Coronaviridae", "Flaviviridae", "Orthomyxoviridae",
                "Retroviridae", "Herpesviridae"]
    mol_types = ["ssRNA(+)", "ssRNA(-)", "dsDNA", "ssRNA", "dsRNA"]

    rows = []
    for host, params in host_categories.items():
        for i in range(n_per_host):
            size = max(1000, int(np.random.normal(
                params["size_mean"], params["size_std"])))
            gc   = round(np.clip(np.random.normal(
                params["gc_mean"], params["gc_std"]), 0.3, 0.7), 3)

            rows.append({
                "Accession":    f"ACC_{len(rows):05d}",
                "Organism Name": f"Virus_{len(rows):04d}",
                "Family":        np.random.choice(families),
                "Host":          host,
                "Genome Length": size,
                "GC Content":    gc,
                "Molecule Type": np.random.choice(mol_types),
                "Release Date":  f"202{np.random.randint(0,5)}-{np.random.randint(1,13):02d}-01"
            })

    # Add SARS-CoV-2 as a known reference
    rows.append({
        "Accession":     "NC_045512",
        "Organism Name": "Severe acute respiratory syndrome coronavirus 2",
        "Family":        "Coronaviridae",
        "Host":          "Human",
        "Genome Length": 29903,
        "GC Content":    0.380,
        "Molecule Type": "ssRNA(+)",
        "Release Date":  "2020-01-17"
    })

    df = pd.DataFrame(rows)
    return df


# ─────────────────────────────────────────────
# 2. LOAD & CLEAN METADATA
# ─────────────────────────────────────────────

def load_metadata(csv_path):
    """
    Load NCBI Virus metadata CSV and standardize column names.

    Expected columns (from NCBI Virus portal download):
        Accession, Organism Name, Family, Host,
        Genome Length, GC Content, Molecule Type

    Parameters
    ----------
    csv_path : str

    Returns
    -------
    df : pd.DataFrame (cleaned)
    """
    df = pd.read_csv(csv_path)
    print(f"[INFO] Loaded {len(df):,} viral records from {csv_path}")
    print(f"[INFO] Columns: {list(df.columns)}")

    # Standardize key column names if needed
    rename_map = {
        "genome_length": "Genome Length",
        "gc_content":    "GC Content",
        "host":          "Host",
        "family":        "Family",
        "molecule_type": "Molecule Type",
    }
    df.rename(columns={k: v for k, v in rename_map.items()
                        if k in df.columns}, inplace=True)

    # Drop rows missing key fields
    df = df.dropna(subset=["Genome Length", "Host"])
    df["Genome Length"] = pd.to_numeric(df["Genome Length"], errors="coerce")
    df = df.dropna(subset=["Genome Length"])

    print(f"[INFO] After cleaning: {len(df):,} records")
    return df


# ─────────────────────────────────────────────
# 3. PLOTS
# ─────────────────────────────────────────────

def plot_genome_size_by_host(df, output_dir):
    """
    Histogram of genome sizes grouped by host category.
    This is the main comparative visualization described in the README.
    """
    hosts = df["Host"].value_counts().head(6).index.tolist()
    df_filtered = df[df["Host"].isin(hosts)]

    palette = sns.color_palette("Set2", len(hosts))

    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharey=False)
    axes = axes.flatten()

    for idx, (host, color) in enumerate(zip(hosts, palette)):
        subset = df_filtered[df_filtered["Host"] == host]["Genome Length"]
        ax = axes[idx]
        ax.hist(subset, bins=20, color=color, edgecolor='white', alpha=0.85)

        # Highlight SARS-CoV-2 genome length if host is Human
        if host == "Human":
            ax.axvline(x=29903, color='red', linestyle='--', linewidth=1.5,
                       label='SARS-CoV-2 (29,903 bp)')
            ax.legend(fontsize=8)

        ax.set_title(f"Host: {host}\n(n={len(subset)})",
                     fontsize=11, fontweight='bold')
        ax.set_xlabel("Genome Length (bp)", fontsize=9)
        ax.set_ylabel("Count", fontsize=9)
        ax.xaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, _: f"{int(x/1000)}k" if x >= 1000 else str(int(x)))
        )

    # Hide unused subplots
    for idx in range(len(hosts), len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle("Viral Genome Size Distribution by Host Category\n"
                 "(Red dashed line = SARS-CoV-2 genome length)",
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    out_path = os.path.join(output_dir, "viral_genome_histogram.png")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Genome size histogram saved: {out_path}")


def plot_gc_content_comparison(df, output_dir):
    """
    Boxplot comparing GC content across host categories.
    SARS-CoV-2 has notably low GC content (~38%).
    """
    if "GC Content" not in df.columns:
        print("[WARN] GC Content column not found, skipping GC plot.")
        return

    hosts = df["Host"].value_counts().head(6).index.tolist()
    df_filtered = df[df["Host"].isin(hosts)].copy()
    df_filtered["GC Content (%)"] = df_filtered["GC Content"] * 100

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=df_filtered, x="Host", y="GC Content (%)",
                palette="Set2", ax=ax, order=hosts)
    sns.stripplot(data=df_filtered, x="Host", y="GC Content (%)",
                  color='black', alpha=0.3, size=3, ax=ax, order=hosts)

    # Reference line for SARS-CoV-2 GC content
    ax.axhline(y=38.0, color='red', linestyle='--', linewidth=1.5,
               label='SARS-CoV-2 GC content (~38%)')
    ax.set_title("GC Content by Host Category\n"
                 "(Red dashed = SARS-CoV-2 reference)",
                 fontsize=13, fontweight='bold')
    ax.set_xlabel("Host Category", fontsize=11)
    ax.set_ylabel("GC Content (%)", fontsize=11)
    ax.legend(fontsize=10)
    plt.tight_layout()
    out_path = os.path.join(output_dir, "gc_content_comparison.png")
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[INFO] GC content plot saved: {out_path}")


def plot_family_distribution(df, output_dir):
    """
    Bar chart of virus family counts in the dataset.
    """
    if "Family" not in df.columns:
        return

    family_counts = df["Family"].value_counts().head(10)

    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.bar(family_counts.index, family_counts.values,
                  color=sns.color_palette("Set3", len(family_counts)),
                  edgecolor='white')

    # Highlight Coronaviridae
    for i, (fam, bar) in enumerate(zip(family_counts.index, bars)):
        if fam == "Coronaviridae":
            bar.set_color('#d73027')
            bar.set_label("Coronaviridae (includes SARS-CoV-2)")

    ax.set_title("Virus Family Distribution in Dataset\n"
                 "(Red = Coronaviridae)",
                 fontsize=13, fontweight='bold')
    ax.set_xlabel("Virus Family", fontsize=11)
    ax.set_ylabel("Number of Records", fontsize=11)
    plt.xticks(rotation=35, ha='right', fontsize=9)
    plt.tight_layout()
    out_path = os.path.join(output_dir, "virus_family_distribution.png")
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[INFO] Family distribution plot saved: {out_path}")


def print_summary_table(df, output_dir):
    """
    Print and save a summary statistics table grouped by host.
    """
    summary = df.groupby("Host")["Genome Length"].agg(
        Count="count",
        Mean_bp="mean",
        Median_bp="median",
        Std_bp="std",
        Min_bp="min",
        Max_bp="max"
    ).round(0).astype(int)

    print("\n[INFO] Genome Size Summary by Host:")
    print(summary.to_string())

    summary.to_csv(os.path.join(output_dir, "viral_genome_summary.csv"))
    print(f"[INFO] Summary table saved.")


# ─────────────────────────────────────────────
# 4. MAIN
# ─────────────────────────────────────────────

def main(metadata_path, use_synthetic, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    print("\n" + "="*50)
    print("  COMPARATIVE VIRAL GENOME ANALYSIS")
    print("="*50)

    if use_synthetic:
        print("[INFO] Using synthetic metadata for demonstration.")
        df = generate_synthetic_metadata(n_per_host=20)
        synthetic_path = os.path.join("data", "viral_metadata_synthetic.csv")
        os.makedirs("data", exist_ok=True)
        df.to_csv(synthetic_path, index=False)
        print(f"[INFO] Synthetic metadata saved: {synthetic_path}")
    else:
        df = load_metadata(metadata_path)

    print_summary_table(df, output_dir)
    plot_genome_size_by_host(df, output_dir)
    plot_gc_content_comparison(df, output_dir)
    plot_family_distribution(df, output_dir)

    print(f"\n[DONE] Comparative analysis complete. All plots in: {output_dir}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Comparative viral genome analysis.")
    parser.add_argument("--metadata",      default="data/viral_metadata.csv",
                        help="Path to NCBI Virus metadata CSV")
    parser.add_argument("--use_synthetic", action="store_true",
                        help="Use built-in synthetic data (for testing)")
    parser.add_argument("--output",        default="results/",
                        help="Output directory for plots")
    args = parser.parse_args()

    main(args.metadata, args.use_synthetic, args.output)
