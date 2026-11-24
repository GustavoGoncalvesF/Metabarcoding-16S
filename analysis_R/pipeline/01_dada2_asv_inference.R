#####################################################################################
# Project: 16S rRNA Metabarcoding
# Script: 01_dada2_asv_inference.R
# Autor: Gustavo Goncalves
# Last edit: 07/10/2026
#
# Description:
# Script to processes filtered paired-end reads using the DADA2 pipeline to:
# - Learn sequencing error rates
# - Dereplicate reads
# - Infer Amplicon Sequence Variants (ASVs)
# - Merge paired-end reads
# - Remove chimeric sequences
# - Generate the ASV table
# - Assign taxonomy
#
# Notes:
#   This script assumes that the primers sequences have already been removed and
#   reads have been length-filtered using minimum lenght.
#
#   Remember to load the config.R file for the specific project before running
#   this script. For example:
#       source("projects/coffee/config.R")
#####################################################################################

# Set up packages
library(dada2); packageVersion("dada2")
library(openxlsx); packageVersion("openxlsx")

# Create directories
dir.create(filtered_path, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_path, recursive = TRUE, showWarnings = FALSE)
dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

# Set up logfile
logfile <- file.path(results_path, "dada2_logfile.log")

log_message <- function(msg) {
  cat(
    paste0(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      " | ",
      msg,
      "\n"
    ),
    file = logfile,
    append = TRUE
  )
}

log_message("Starting dada2_asv_inference script from beginning.")

# Forward and reverse FASTQ files
fnFs <- sort(
  list.files(
    raw_path,
    pattern = forward_pattern,
    full.names = TRUE
  )
)

fnRs <- sort(
  list.files(
    raw_path,
    pattern = reverse_pattern,
    full.names = TRUE
  )
)

# Extract sample names
sample.names <- sub(
  forward_pattern,
  "",
  basename(fnFs)
)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(
  filtered_path,
  paste0(sample.names, "_F_filt.fastq.gz")
)

filtRs <- file.path(
  filtered_path,
  paste0(sample.names, "_R_filt.fastq.gz")
)

names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim
log_message("Starting filterAndTrim.")

out <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  maxN = maxN,
  maxEE = maxEE,
  truncQ = truncQ,
  truncLen = truncLen,
  rm.phix = rm_phix,
  minLen = minLen,
  compress = TRUE,
  multithread = n_threads
)

saveRDS(
  out,
  file.path(objects_path, "filterAndTrim_out.rds")
)

log_message("Finished filterAndTrim.")

# Learn error rates
log_message("Starting learnErrors F and R.")

errF <- learnErrors(
  filtFs,
  multithread = n_threads
)

errR <- learnErrors(
  filtRs,
  multithread = n_threads
)

saveRDS(
  errF,
  file.path(objects_path, "errF.rds")
)

saveRDS(
  errR,
  file.path(objects_path, "errR.rds")
)

log_message("Finished learnErrors F and R.")

# Sample inference
log_message("Starting dada F and R.")

dadaFs <- dada(
  filtFs,
  err = errF,
  multithread = n_threads
)

dadaRs <- dada(
  filtRs,
  err = errR,
  multithread = n_threads
)

saveRDS(
  dadaFs,
  file.path(objects_path, "dadaFs.rds")
)

saveRDS(
  dadaRs,
  file.path(objects_path, "dadaRs.rds")
)

log_message("Finished dada F and R.")

# Merge paired reads
log_message("Starting mergePairs.")

mergers <- mergePairs(
  dadaFs,
  filtFs,
  dadaRs,
  filtRs,
  verbose = TRUE
)

log_message("Finishe mergePairs.")

# Construct Amplicon Sequence Variant (ASV) table
log_message("Starting makeSequenceTable.")

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

saveRDS(
  seqtab,
  file.path(objects_path, "seqtab.rds")
)

log_message("Finished makeSequenceTable.")

# Remove chimeras
log_message("Starting removeBimeraDenovo.")

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = n_threads,
  verbose = TRUE
)

dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)

saveRDS(
  seqtab.nochim,
  file.path(objects_path, "seqtab-nochim.rds")
)

log_message("Finished removeBimeraDenovo.")

# Read tracking
source("pipeline/dada2_read_tracking.R")

tracking <- generate_read_tracking(
  out = out,
  dadaFs = dadaFs,
  dadaRs = dadaRs,
  mergers = mergers,
  seqtab.nochim = seqtab.nochim,
  sample.names = sample.names,
  results_path = results_path
)

# Assign taxonomy
log_message("Starting assignTaxonomy.")

taxa <- assignTaxonomy(
  seqtab.nochim,
  taxonomy_database,
  minBoot = min_boot,
  multithread = n_threads
)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL

head(taxa.print)

saveRDS(
  taxa,
  file.path(objects_path, "taxonomy.rds")
)

log_message("Finished assignTaxonomy.")

log_message("End of dada2_asv_inference script.")
