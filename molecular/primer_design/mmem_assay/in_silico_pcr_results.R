
##############################################
# In-silico PCR for Membranipora membranacea
# Author: Calum + ChatGPT
# Requires: Biostrings, tidyverse
##############################################

# Install required packages if needed
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Biostrings")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")

library(Biostrings)
library(tidyverse)

# -----------------------------
# 1. Input
# -----------------------------
# Path to genome fasta
genome_file <- "mmem_genome.fna"

# Primers
fwd_primer <- "ACGGATGAGTTGTTAGCCCA"
rev_primer <- "GAATGGAAGCCTGCAACACA"   # given in 5'→3'

# Parameters
max_mismatch <- 2        # allow up to 2 mismatches
min_size <- 100          # minimum expected product size
max_size <- 1000         # maximum expected product size

# -----------------------------
# 2. Read genome
# -----------------------------
cat("Loading genome...\n")
genome <- readDNAStringSet(genome_file)
cat(length(genome), "contigs loaded\n")

# -----------------------------
# 3. Search for primers
# -----------------------------
cat("Searching for primer matches...\n")

# Forward primer matches
fwd_hits <- vmatchPattern(fwd_primer, genome, max.mismatch = max_mismatch)

# Reverse complement primer matches
rev_rc <- as.character(reverseComplement(DNAString(rev_primer)))
rev_hits <- vmatchPattern(rev_rc, genome, max.mismatch = max_mismatch)

# -----------------------------
# 4. Extract hits and pair them
# -----------------------------
get_hits_df <- function(hits, strand_label) {
  tibble(
    contig = names(hits),
    hits = lapply(hits, as.data.frame)
  ) %>%
    unnest(hits) %>%
    mutate(strand = strand_label) %>%
    select(contig, start, end, width, strand)
}

fwd_df <- get_hits_df(fwd_hits, "+")
rev_df <- get_hits_df(rev_hits, "-")

# -----------------------------
# 5. Identify matching forward–reverse pairs
# -----------------------------
cat("Pairing primer hits...\n")

pairs <- inner_join(fwd_df, rev_df, by = "contig", suffix = c("_fwd", "_rev")) %>%
  mutate(product_size = start_rev - end_fwd + 1) %>%
  filter(product_size >= min_size, product_size <= max_size) %>%
  arrange(product_size)

# -----------------------------
# 6. Output
# -----------------------------
if (nrow(pairs) == 0) {
  cat("⚠️  No valid amplicons found between", min_size, "and", max_size, "bp.\n")
} else {
  cat("✅ Found", nrow(pairs), "potential amplicons:\n")
  print(pairs %>%
    select(contig, start_fwd, end_rev, product_size) %>%
    head(10))
  
  write.csv(pairs, "in_silico_pcr_results.csv", row.names = FALSE)
  cat("\nResults saved to 'in_silico_pcr_results.csv'\n")
}




##############################################
# Visualize in-silico PCR hits
##############################################
library(tidyverse)
library(ggplot2)

# Read results
hits <- read.csv("in_silico_pcr_results.csv")

# Plot positions
ggplot(hits, aes(y = contig)) +
  geom_segment(aes(x = start_fwd, xend = end_rev), color = "dodgerblue3", size = 2) +
  geom_point(aes(x = start_fwd), color = "forestgreen", size = 3) +
  geom_point(aes(x = end_rev), color = "firebrick3", size = 3) +
  theme_minimal() +
  labs(
    title = "In-silico PCR Products in Membranipora membranacea Genome",
    x = "Genomic position (bp)",
    y = "Contig",
    caption = "Blue = predicted amplicon | Green = forward primer | Red = reverse primer"
  )



hits <- read.csv("in_silico_pcr_results.csv")
hits %>%
  count(product_size) %>%
  arrange(product_size)



#########################################
# Extract in-silico PCR amplicons from genome
# Requires: Biostrings, tidyverse
#########################################

# Install if needed
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("Biostrings")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")

library(Biostrings)
library(tidyverse)

# 1. Load data
genome_file <- "mmem_genome.fna"
hits_file   <- "in_silico_pcr_results.csv"

genome <- readDNAStringSet(genome_file)
hits   <- read.csv(hits_file)

# 2. Filter for short amplicons (eDNA-suitable)
target_hits <- hits %>%
  filter(product_size %in% c(129, 159))   # or modify to your preferred lengths

# 3. Extract each region from the genome
extract_amplicon <- function(contig, start, end) {
  seq <- genome[[contig]]
  subseq(seq, start = start, end = end)
}


############################################
# Fix contig naming and safely extract amplicons
############################################

# Clean contig names: keep only accession up to first space
hits$contig <- sub(" .*", "", hits$contig)

# Check what contig names exist in genome
cat("Example genome contig names:\n")
print(head(names(genome)))

# Filter for your target product sizes
target_hits <- hits %>%
  filter(product_size %in% c(129, 159))

# Extract safely with conditional check
amplicons <- target_hits %>%
  rowwise() %>%
  mutate(sequence = if (contig %in% names(genome)) {
    as.character(subseq(genome[[contig]], start_fwd, end_rev))
  } else {
    NA_character_   # returns NA if contig not found
  }) %>%
  ungroup() %>%
  filter(!is.na(sequence))

# 5. Write FASTA output
fasta_lines <- c()
for (i in seq_len(nrow(amplicons))) {
  fasta_lines <- c(
    fasta_lines,
    paste0(">amplicon_", i, "_", amplicons$contig[i], "_", amplicons$product_size[i], "bp"),
    amplicons$sequence[i]
  )
}
writeLines(fasta_lines, "mmem_candidate_amplicons_fixed.fasta")

cat("✅ Amplicon sequences saved to 'mmem_candidate_amplicons_fixed.fasta'\n")
