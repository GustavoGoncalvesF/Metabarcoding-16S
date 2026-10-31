#####################################################################################
# Project: Paspalum Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: config.R
#
# Description:
# Project-specific paths and parameters for the DADA2 pipeline.
#
#####################################################################################

# Project information
project_name <- "paspalum_metabarcoding"

project_path <- "projects/paspalum"

raw_path <- file.path(project_path, "raw")
filtered_path <- file.path(project_path, "filtered")
results_path <- file.path(project_path, "results")
objects_path <- file.path(results_path, "objects")
tax_path <- file.path(project_path, "tax")

# FASTQ filename patterns
forward_pattern <- "_R1.final.fastq.gz$"
reverse_pattern <- "_R2.final.fastq.gz$"

# DADA2 filtering
maxN     <- 0
maxEE    <- c(2, 2)
truncQ   <- 2
minLen   <- 150
rm_phix  <- TRUE
truncLen <- c(0, 0)

# Computational resources
n_threads <- 18  # On Mac Pro Mid, set max 24 threads

# Taxonomy
taxonomy_database <- file.path(
  tax_path,
  "silva_nr99_v138.2_toGenus_trainset.fa.gz"
)
min_boot <- 50

# DADA2 QC
expected_min_asv_length <- 400

