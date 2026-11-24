#####################################################################################
# Project: Paspalum Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 03_phyloseq_object.R
# Author: Gustavo Goncalves
# Last edit: 09/11/2026
#
# Description:
# Script to create and adjust the phyloseq object for the paspalum rhizosphere
# microbiome project.
#
#####################################################################################

# Load project configuration
source("projects/paspalum/config.R")

# Set up packages
library(phyloseq); packageVersion("phyloseq")

# --------------------------------------------------
# Load DADA2 outputs
# --------------------------------------------------
seqtab.nochim <- readRDS(
  file.path(objects_path, "seqtab-nochim.rds")
)

taxa <- readRDS(
  file.path(objects_path, "taxonomy.rds")
)

# Load sample metadata
metadata <- readxl::read_excel(
  file.path(project_path, "metadata", "sample_data.xlsx")
)

# Convert tibble to data.frame
metadata <- as.data.frame(metadata)

# Standardize SampleID
metadata$SampleID <- trimws(as.character(metadata$SampleID))

# Convert P1B1.1 to P1B1-1
metadata$SampleID <- sub("\\.([0-9]+)$", "-\\1", metadata$SampleID)

# Assign row names
rownames(metadata) <- metadata$SampleID

# --------------------------------------------------
# Validate sample and ASV identifiers
# --------------------------------------------------

required_columns <- c(
  "SampleID", "Replicate", "Treatment", "Species", "Accession", "Origin"
)

missing_columns <- setdiff(required_columns, colnames(metadata))
if (length(missing_columns) > 0) {
  stop("Missing metadata columns: ", paste(missing_columns, collapse = ", "))
}

# Normalize sample IDs to match DADA2: P1B1.1 -> P1B1-1
metadata$SampleID <- trimws(as.character(metadata$SampleID))
metadata$SampleID <- sub("\\.([0-9]+)$", "-\\1", metadata$SampleID)

if (anyNA(metadata$SampleID) || any(metadata$SampleID == "") ||
    anyDuplicated(metadata$SampleID)) {
  stop("SampleID contains missing, empty or duplicated identifiers.")
}

rownames(metadata) <- metadata$SampleID

if (anyDuplicated(rownames(seqtab.nochim))) {
  stop("Duplicated sample IDs in the ASV table.")
}

if (!setequal(rownames(metadata), rownames(seqtab.nochim))) {
  stop(
    "Sample IDs in metadata do not match ASV table.\n",
    "Only in metadata: ",
    paste(setdiff(rownames(metadata), rownames(seqtab.nochim)), collapse = ", "),
    "\nOnly in ASV table: ",
    paste(setdiff(rownames(seqtab.nochim), rownames(metadata)), collapse = ", ")
  )
}

# Reorder metadata to match DADA2 sample order
metadata <- metadata[rownames(seqtab.nochim), , drop = FALSE]

if (anyDuplicated(colnames(seqtab.nochim)) || anyDuplicated(rownames(taxa))) {
  stop("Duplicated ASV identifiers in the ASV or taxonomy table.")
}

if (!setequal(colnames(seqtab.nochim), rownames(taxa))) {
  stop("ASV sequences in the ASV table do not match the taxonomy table.")
}

# Reorder taxonomy to match DADA2 ASV order
taxa <- taxa[colnames(seqtab.nochim), , drop = FALSE]

# --------------------------------------------------
# Prepare sample metadata
# --------------------------------------------------

metadata$Replicate <- as.integer(metadata$Replicate)
if (anyNA(metadata$Replicate)) {
  stop("Replicate contains missing or invalid integer values.")
}

metadata$Treatment <- factor(metadata$Treatment)
metadata$Species <- factor(metadata$Species)
metadata$Accession <- factor(metadata$Accession)
metadata$Origin <- factor(metadata$Origin, levels = c("BAG1", "BAG2"))

if (anyNA(metadata$Origin)) {
  stop("Origin must contain BAG1 or BAG2 for every sample.")
}

# --------------------------------------------------
# Validate metadata
# --------------------------------------------------

stopifnot(
  identical(rownames(metadata), rownames(seqtab.nochim)),
  identical(colnames(seqtab.nochim), rownames(taxa))
)

cat("\nSamples per treatment:\n")
print(table(metadata$Treatment, useNA = "ifany"))

cat("\nSamples per species and origin:\n")
print(table(metadata$Species, metadata$Origin, useNA = "ifany"))

cat("\nSamples per accession and origin:\n")
print(table(metadata$Accession, metadata$Origin, useNA = "ifany"))

cat("\nMissing metadata values:\n")
print(colSums(is.na(metadata)))

# --------------------------------------------------
# Construct complete phyloseq object
# --------------------------------------------------

ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(metadata)
)

# Preserve original nucleotide sequences in refseq before renaming ASVs
asv_sequences <- taxa_names(ps)
dna <- Biostrings::DNAStringSet(asv_sequences)
names(dna) <- asv_sequences
ps <- merge_phyloseq(ps, dna)

# Rename ASVs consistently across otu_table, tax_table, and refseq
taxa_names(ps) <- paste0("ASV", seq_len(ntaxa(ps)))

saveRDS(
  ps,
  file.path(objects_path, "phyloseq_object.rds")
)
