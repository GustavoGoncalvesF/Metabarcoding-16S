#####################################################################################
# Project: Coffe Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: phyloseq_analysis.R
# Autor: Gustavo Goncalves
# Last edit: 21/07/2026
#
# Description:
# Script to create a phyloseq object and run some statistical analysis:
# - Visualize alpha-diversity
#
# Notes:
#    Phyloseq object can be loaded by load_objects.R)
#####################################################################################

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(vegan); packageVersion("vegan")

# Set up paths
path <- "results"; list.files(path)

# Set up logfile
logfile <- "phyloseq_logfile.log"
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
log_message("Starting phyloseq_analysis script from beginning.")

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
log_message("Finished phyloseq object.")

theme_set(theme_bw()) # Set ggplot2 theme

# ==================================================================================
# Exploratory analysis of alpha-diversity
# ==================================================================================
# Observed + Shannon - Farm x CropStage
plot_richness(
  ps,
  x = "Farm",
  measures = c("Shannon", "Observed"),
  color = "CropStage"
)

# (Observed + Shannon - Farm x CropStage) + facet Year
plot_richness(
  ps,
  x = "Farm",
  measures = c("Shannon", "Observed"),
  color = "CropStage"
) +
  facet_wrap(~Year)

# Observed + Shannon - Farm x pH
plot_richness(
  ps,
  x="Farm",
  measures="Shannon",
  color="pH"
)

# ==================================================================================
# Relative abundance between samples
# ==================================================================================
# Obtain relative abundance
ps.rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)

# Top phylum
ps.phylum <- tax_glom(ps.rel, taxrank = "Phylum")

top10 <- names(
  sort(
    taxa_sums(ps.phylum),
    decreasing = TRUE
  )[1:10]
)

ps.phylum.top10 <- prune_taxa(
  top10,
  ps.phylum
)

# SampleID
plot_bar(
  ps.phylum.top10,
  x = "SampleID",
  fill = "Phylum"
)

# SampleID per Farm
plot_bar(
  ps.phylum.top10,
  x = "SampleID",
  fill = "Phylum"
) +
  facet_wrap(~Farm, scales = "free_x")

# Facet Farms
plot_bar(
  ps.phylum.top10,
  x = "CropStage",
  fill = "Phylum"
) +
  facet_wrap(~Farm)

# ==================================================================================
# Exploratory analysis of beta-diversity
# ==================================================================================

# Bray-Curtis (distance matrix)
bray <- phyloseq::distance(
  physeq = ps,
  method = "bray"
)

# PCoA ordinate
ord.pcoa <- ordinate(
  ps,
  method = "PCoA",
  distance = "bray"
)

# Plot by CropStage
plot_ordination(
  ps,
  ord.pcoa,
  color = "CropStage"
) +
  geom_point(size = 3)

# Plot by Quality
plot_ordination(
  ps,
  ord.pcoa,
  color = "Quality"
) +
  geom_point(size = 3)

# Plot by Farm
plot_ordination(
  ps,
  ord.pcoa,
  color = "Farm"
) +
  geom_point(size = 3)

# Plot by Farm + groups
plot_ordination(
  ps,
  ord.pcoa,
  color = "Farm"
) +
  geom_point(size = 3) +
  stat_ellipse()

# ==================================================================================
# PERMANOVA
# ==================================================================================
adonis2(bray ~ Farm, data = data.frame(sample_data(ps))) # Farm factor
adonis2(bray ~ CropStage, data = data.frame(sample_data(ps))) # CropStage factor
adonis2(bray ~ Quality, data = data.frame(sample_data(ps))) #Quality factor
adonis2(bray ~ Farm + CropStage + Quality, data = data.frame(sample_data(ps)), by = "margin")

# beta-dispersion test
disp.farm <- betadisper(bray, sample_data(ps)$Farm)
permutest(disp.farm)
boxplot(disp.farm) # Visualize beta-dispersion
TukeyHSD(disp.farm) # Show differences between farms

# ==================================================================================
# Environmental analysis
# ==================================================================================

# Environment data frame
env <- data.frame(sample_data(ps))[
  ,
  c(
    "pH",
    "Phosphorus",
    "Potassium",
    "CEC",
    "Altitude"
  )
]
pcoa_scores <- as.data.frame(ord.pcoa$vectors) # Extract scores from ord.pcoa

# envfit to environmental factors analysis
ef <- envfit(
  pcoa_scores[,1:2],
  env,
  permutations = 999
); ef

log_message("End of phyloseq_analysis script.")
