#####################################################################################
# Project: Coffe Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: dada2_asv_inference.R
# Autor: Gustavo Goncalves
# Last edit: 16/07/2026
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
#####################################################################################

# Set up packages
library(dada2);packageVersion("dada2")
library(stats);packageVersion("stats")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

# Set up paths
path <- "raw"; list.files(path)
filtpath <- "filtered"
resultspath <- "results"

# Set up logfile
logfile <- "dada2_logfile.log"
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

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.final.fastq.gz and SAMPLENAME_R2.final.fastq.gz
fnFs <- sort(list.files(path, pattern="_R1.final.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.final.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(filtpath, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtpath, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim
log_message("Starting filterAndTrim.")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, minLen=(150),
                     compress=TRUE, multithread=18) # On Mac Pro Mid, set max 24 threads
saveRDS(out, file.path(resultspath, "filterAndTrim_out.rds"))
log_message("Finished filterAndTrim.")

# Learn error rates
log_message("Starting learnErrors F and R.")
errF <- learnErrors(filtFs, multithread=18)
errR <- learnErrors(filtRs, multithread=18)
saveRDS(errF, file.path(resultspath, "errF.rds")) # Save error rates forward
saveRDS(errR, file.path(resultspath, "errR.rds")) # Save error rates reverse
log_message("Finished learnErrors F and R.")

# Sample inference
log_message("Starting dada F and R.")
dadaFs <- dada(filtFs, err=errF, multithread=18)
dadaRs <- dada(filtRs, err=errR, multithread=18)
saveRDS(dadaFs, file.path(resultspath, "dadaFs.rds")) # Save dada file forward
saveRDS(dadaRs, file.path(resultspath, "dadaRs.rds")) # Save dada file reverse
log_message("Finished dada F and R.")

# Merge paired reads
log_message("Starting mergePairs.")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
log_message("Finishe mergePairs.")

# Construct Amplicon Sequence Variant (ASV) table
log_message("Starting makeSequenceTable.")
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
saveRDS(seqtab, file.path(resultspath, "seqtab.rds")) # Save ASV table
log_message("Finished makeSequenceTable.")

# Remove chimeras
log_message("Starting removeBimeraDenovo.")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=18, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim, file.path(resultspath, "seqtab-nochim.rds")) # Save ASV table without chimeras
log_message("Finished removeBimeraDenovo.")

# ASVs Summary
summary(colSums(seqtab.nochim)) # Abundance distribution of ASVs
table(colSums(seqtab.nochim) == 1) # Verify singletons
apply(seqtab.nochim > 0, 1, sum) # Richness per sample

# Track reads throught the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, file.path(resultspath, "read_tracking.csv")) # Save track table

# Assign taxonomy
log_message("Starting assignTaxonomy.")
taxa <- assignTaxonomy(seqtab.nochim, "tax/silva_nr99_v138.2_toGenus_trainset.fa.gz", minBoot = 50, multithread=18)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
saveRDS(taxa, file.path(resultspath, "taxonomy.rds")) # Save taxonomy object
log_message("Finished assignTaxonomy.")

# Taxonomy summary
colSums(!is.na(taxa)) # ASVs assigned in each level
round(100 * colSums(!is.na(taxa)) / nrow(taxa),2)

head(sort(table(taxa[, "Phylum"]), decreasing = TRUE), 15) # Main phylums
head(sort(table(taxa[, "Class"]), decreasing=TRUE), 15) # Main classes
head(sort(table(taxa[, "Genus"]), decreasing=TRUE), 20) # Main genus

theme_set(theme_bw())

# Import metadata from sample_data
metadata <- read.csv(
  "metadata/sample_data.csv",
  sep = ";",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)
all(rownames(metadata) == rownames(seqtab.nochim)) # Verify sample names

# Convert variable types
metadata$Replicate <- factor(metadata$Replicate)
metadata$Treatment <- factor(metadata$Treatment)
metadata$Farm <- factor(metadata$Farm)
metadata$CoffeeSpecies <- factor(metadata$CoffeeSpecies)
metadata$Cultivar <- factor(metadata$Cultivar)
metadata$Quality <- factor(metadata$Quality)
metadata$CropStage <- factor(metadata$CropStage)
metadata$SampleType <- factor(metadata$SampleType)
metadata$Year <- factor(metadata$Year)
metadata$Altitude <- as.numeric(metadata$Altitude)
metadata$Potassium <- as.numeric(metadata$Potassium)
metadata$Phosphorus <- as.numeric(metadata$Phosphorus)
metadata$pH <- as.numeric(metadata$pH)
metadata$CEC <- as.numeric(metadata$CEC)

# Construct phyloseq object
log_message("Making phyloseq object.")
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(metadata)
)
saveRDS(ps, file.path(resultspath, "phyloseq_object.rds")) # Save phyloseq object
log_message("Finished phyloseq object.")

log_message("End of dada2_asv_inference script.")
