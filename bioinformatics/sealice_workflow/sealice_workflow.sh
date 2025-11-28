#!/bin/bash

# ===== STEP 1: SIMULATE NANOPORE READS =====
# Activate Conda env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate badread_env

# Simulate reads
badread simulate \
  --reference GCF_016086655.3_UVic_Lsal_1.2_genomic.fna \
  --quantity 50x \
  --length 10000,5000 \
  --error_model nanopore \
  --output nanopore_reads.fastq

# ===== STEP 2: FIND REPEATS =====
# Make RepeatModeler DB
BuildDatabase -name mydb GCF_016086655.3_UVic_Lsal_1.2_genomic.fna

# Run RepeatModeler
RepeatModeler -database mydb -pa 4

# Mask repeats
RepeatMasker -lib mydb-families.fa my_genome.fna

# ===== STEP 3: BLAST TO CHECK EUROPEAN ORIGIN =====
# Make BLAST database (replace with your actual reference)
makeblastdb -in european_reference.fasta -dbtype nucl

# Run BLAST
blastn -query GCF_016086655.3_UVic_Lsal_1.2_genomic.fna -db european_reference.fasta -out blast_results.txt -outfmt 6

echo "✅ Workflow complete!"
