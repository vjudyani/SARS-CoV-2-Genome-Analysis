# SARS-CoV-2 Genome Analysis

A bioinformatics pipeline for analyzing the SARS-CoV-2 genome — identifying protein-coding genes (ORFs), predicting membrane-associated proteins via hydrophobicity analysis, and comparing SARS-CoV-2 to other viruses across host categories.

---

## Project Overview

SARS-CoV-2, the virus responsible for COVID-19, is studied using bioinformatics to understand how it infects cells, replicates, and evades the immune system. This project:

1. **Downloads** the SARS-CoV-2 genome (NC_045512) directly from NCBI
2. **Validates** the ORF-finding algorithm on E. coli as a well-annotated ground truth
3. **Identifies ORFs** in all 6 reading frames of the SARS-CoV-2 genome
4. **Analyzes protein hydrophobicity** using the Kyte-Doolittle scale to predict membrane proteins and drug targets
5. **Compares** SARS-CoV-2 genome properties against other viruses across host categories

---

## Repository Structure

```
SARS-CoV-2-Genome-Analysis/
│
├── data/                                  # Input data (not tracked by git)
│   ├── sars_cov2.fasta                    # Downloaded from NCBI (NC_045512)
│   ├── ecoli.fasta                        # E. coli K-12 genome (U00096)
│   └── viral_metadata.csv                 # NCBI Virus metadata (see instructions)
│
├── results/                               # Auto-generated outputs
│   ├── sars_cov2_orfs.csv                 # ORF table (position, strand, length)
│   ├── sars_cov2_proteome.fasta           # Predicted protein sequences
│   ├── sars_cov2_hydrophobicity.csv       # Hydrophobicity scores per ORF
│   ├── sars_cov2_hydrophobicity_distribution.png
│   ├── sars_cov2_top_hydrophobic_orfs.png
│   ├── sars_cov2_length_vs_hydrophobicity.png
│   ├── viral_genome_histogram.png
│   ├── gc_content_comparison.png
│   └── virus_family_distribution.png
│
├── 01_download_genomes.py                 # Download SARS-CoV-2 and E. coli from NCBI
├── 02_orf_prediction.py                   # ORF detection in all 6 reading frames
├── 03_hydrophobicity.py                   # Kyte-Doolittle hydrophobicity analysis
├── 04_comparative_analysis.py             # Comparative viral genome analysis
├── requirements.txt
└── README.md
```

---

## Background

### Why E. coli for Validation?

E. coli K-12 is one of the most well-annotated organisms in existence, with around 4,300 confirmed protein-coding genes. Before applying the ORF-finding algorithm to SARS-CoV-2, it is first tested on E. coli — since the ground truth is known, we can measure false positive rates and calibrate the minimum ORF length threshold accordingly.

### Why Hydrophobicity?

Proteins with high hydrophobicity scores tend to be membrane-associated, including the SARS-CoV-2 **Spike protein** (target of COVID-19 vaccines), **Envelope protein**, and **Membrane protein**. Identifying these from sequence alone provides candidate drug and vaccine targets without requiring expensive experimental assays.

---

## Pipeline Steps

### Step 1 — Download Genomes

```bash
python 01_download_genomes.py --email your.email@example.com
```

Downloads SARS-CoV-2 (NC_045512) and E. coli (U00096) from NCBI Entrez and saves them as FASTA files in the `data/` directory.

### Step 2 — ORF Prediction

```bash
python 02_orf_prediction.py \
    --genome data/sars_cov2.fasta \
    --organism sars_cov2 \
    --min_length 100 \
    --validate
```

Scans all 6 reading frames (3 forward, 3 reverse complement) and filters out ORFs shorter than 100 nucleotides (33 amino acids). The `--validate` flag runs E. coli first as a sanity check. Outputs an ORF table CSV and a protein FASTA file.

### Step 3 — Hydrophobicity Analysis

```bash
python 03_hydrophobicity.py \
    --fasta results/sars_cov2_proteome.fasta \
    --organism sars_cov2
```

Computes mean Kyte-Doolittle hydrophobicity per protein. Generates distribution plots and identifies the top membrane protein candidates.

### Step 4 — Comparative Viral Analysis

```bash
# With your own NCBI Virus metadata:
python 04_comparative_analysis.py --metadata data/viral_metadata.csv

# Or use the built-in synthetic data to test the pipeline:
python 04_comparative_analysis.py --use_synthetic
```

---

## Downloading Viral Metadata (Step 4)

1. Go to [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/)
2. Search for viruses and select multiple host categories
3. Download a randomized subset (around 20 records per host category) as CSV
4. Rename the file to `viral_metadata.csv` and place it in the `data/` folder

---

## Key Results

### ORF Prediction — SARS-CoV-2

| Metric | Value |
|--------|-------|
| Genome length | 29,903 bp |
| Total ORFs identified (>=100 nt) | ~40-60 candidates |
| Known annotated proteins | 29 (Spike, N, M, E, nsp1-16, ORF3a, etc.) |

### Hydrophobicity — Top Candidates

Highly hydrophobic ORFs (score > 1.0) are candidates for membrane-associated proteins including the **Spike protein (S)**, which mediates ACE2 receptor binding and membrane fusion; the **Envelope protein (E)**, a viroporin and vaccine/drug target; and the **Membrane protein (M)**, the most abundant structural protein in the virion.

### Comparative Analysis

SARS-CoV-2 (Coronaviridae) has one of the largest RNA virus genomes at roughly 30 kb — about three times the size of influenza. Its low GC content (~38%) is notably lower than most bat coronaviruses (~40-45% GC), which is one of the genomic features explored in the comparative analysis.

---

## Requirements

```
biopython>=1.79
pandas>=1.5.0
numpy>=1.23.0
matplotlib>=3.5.0
seaborn>=0.12.0
```

Install all dependencies:

```bash
pip install -r requirements.txt
```

---

## Future Work

- Add BLAST homology search against known viral proteins
- Integrate multiple sequence alignment (MUSCLE/MAFFT) for variant comparison
- Predict transmembrane helices using TMHMM or DeepTMHMM
- Add phylogenetic tree construction comparing SARS-CoV-2 to related betacoronaviruses
- Automate the full pipeline with Snakemake for reproducibility

---

## Author

**Vedika Judyani**  
MS Bioinformatics | Bioinformatics Analyst  
[LinkedIn](https://www.linkedin.com/in/vedika-judyani-a19011128/) | [GitHub](https://github.com/vjudyani)

---

## License

This project is licensed under the [CC BY-NC-SA 4.0 License](https://creativecommons.org/licenses/by-nc-sa/4.0/).

---

## References

- Wu et al. (2020). A new coronavirus associated with human respiratory disease in China. *Nature*, 579, 265-269.
- Kyte & Doolittle (1982). A simple method for displaying the hydropathic character of a protein. *J Mol Biol*, 157(1), 105-132.
- [NCBI Virus Portal](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/)
- [Biopython Documentation](https://biopython.org/docs/latest/api/)
