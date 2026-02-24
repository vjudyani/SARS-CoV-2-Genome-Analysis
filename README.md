SARS-CoV-2 Genome Analysis

## Introduction

SARS-CoV-2, the virus responsible for COVID-19, is being studied using bioinformatics to better understand how it functions and to find potential treatment targets. We are focusing on identifying the genes in its genome that help it infect cells, replicate, and evade the immune system. By comparing these genes to those of other viruses, we can predict their roles and identify possible drug targets.

To do this, we look for sections of the genome that might code for proteins, known as open reading frames (ORFs). We filter out any that are too short to be functional and verify our predictions using specialized algorithms. Before applying this method to SARS-CoV-2, we first test it on *E. coli*, a well-studied bacterium, to make sure our approach is accurate. This step-by-step process ensures reliable results while reducing errors that would need further experimental testing.

This project aims to provide a structured bioinformatics pipeline for analyzing the SARS-CoV-2 genome, identifying key protein-coding genes, and exploring their functional roles.

Escherichia coli (E. Coli) is a bacteria commonly found in the human intestine. Most strains are harmless to humans, at worse, causing food poisoning and diarrhea. They can survive outside a host for only a short period of time, making it a potential indicator of fecal contamination. Over the years, E. Coli has been intensely studied and is probably one of the most well-understood organisms in existence. We've learned how to grow them in an optimal environment where they can reproduce up to once every 20 minutes. Due to their rapid growth and easy manipulation, biologists often use them to produce recombinant proteins.

Recombinant proteins are proteins that wouldn't naturally appear in that organism. For instance, we can insert genes that code for fluorescence into plants, making them glow in the dark. Or, perhaps more usefully, we can take the human gene that codes for insulin and convince E. Coli to produce insulin instead. Insulin that we can then use to treat diabetic patients. Or proteins used in cancer treatment. Or, more recently, we can insert fragments of the SARS-CoV-2 virus into E. Coli and use that to produce COVID-19 vaccines. E. coli is a very well-studied organism, therefore, well-annotated. This means we can quickly check our work for any analysis we might perform because we have the ground truth, which biologists have spent decades meticulously gathering for us. We will examine the DNA sequence of E. coli and implement an algorithm for finding potential gene candidates. Because the ground truth is readily available, we can check how many of our ORF candidates are actual genes and how many candidates are false positives.


Setup Instructions

1. Set Up Environment

Create a working directory:

mkdir sars_cov2_analysis
cd sars_cov2_analysis

2. Install Required Packages

Ensure Python and necessary dependencies are available:

module load python 
pip install biopython matplotlib jupyterlab seaborn pandas

Gene Prediction and Validation

Download and Analyze the SARS-CoV-2 Genome

Fetch the SARS-CoV-2 genome using Biopython:

from Bio import Entrez
Entrez.email = "your.email@charlotte.edu"
handle = Entrez.efetch(db="nucleotide", id="NC_045512", rettype="fasta", retmode="text")
with open("sars_cov2.fasta", "w") as output_file:
    output_file.write(handle.read())
handle.close()

Identify Open Reading Frames (ORFs)

Use the following function to identify ORFs in the SARS-CoV-2 genome:

from Bio import SeqIO

def find_orfs(seq_record, min_length=30):
    """Identify ORFs in a given sequence."""
    ... # ORF detection code here

Protein Hydrophobicity Analysis

Generate a Protein FASTA File

Export ORFs as protein sequences:

def export_trimmed_proteins(orfs, output_filename):
    """Save protein sequences in FASTA format."""
    ... # Protein export function here

Compute Hydrophobicity

Load hydrophobicity values and analyze protein sequences:

import pandas as pd
hydro_table = load_hydrophobicity_table("aminoacid_properties.csv")
sars_cov2_data = calculate_hydrophobicities("sars_cov2_proteome.fasta", hydro_table)
df = pd.DataFrame(sars_cov2_data)

Comparative Viral Genome Analysis

Download Viral Metadata

Visit NCBI Virus

Download a randomized subset (20 records per host category) as a CSV file

Rename it to viral_metadata.csv and upload it to the working directory

Generate a Comparative Plot

python viral_genome_histogram.py

License

This project is licensed under the CC BY-NC-SA 4.0 License.
