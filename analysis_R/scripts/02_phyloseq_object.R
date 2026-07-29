#####################################################################################
# Project: Coffe Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 02_phyloseq_object.R
# Autor: Gustavo Goncalves
# Last edit: 29/07/2026
#
# Description:
# Script to create a phyloseq object
#
#####################################################################################

# Set up packages
library(phyloseq); packageVersion("phyloseq")

# Set up paths
path <- "results"; list.files(path)

# Load ASVs and taxonomy matrix
seqtab.nochim <- readRDS(file.path(path, "seqtab-nochim.rds"))
taxa <- readRDS(file.path(path, "taxonomy.rds"))

# Load metadata from sample_data
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

# Adjust
sample_data(ps)$SampleID <- sample_names(ps) # Add SampleID column

# Rename taxa to short names of ASVs and store the DNA sequences in refseq slot
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

saveRDS(ps, file.path(path, "phyloseq_object.rds")) # Save phyloseq object
