#~~~~~~~ ASVs table construction and taxonomy assign by DADA2 ~~~~~~~#
library(dada2);packageVersion("dada2")
library(stats);packageVersion("stats")

path <- "raw"  # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
filtpath <- "filtered"
resultspath <- "results"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.final.fastq.gz and SAMPLENAME_R2.final.fastq.gz
fnFs <- sort(list.files(path, pattern="_R1.final.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.final.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles of the forward reads
plotQualityProfile(fnFs[1:2])

# Inspect read quality profiles of the reverse reads
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(filtpath, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtpath, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, minLen=(150),
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE (only needed for filterAndTrim)
head(out)

# Save and load out
saveRDS(out, file.path(resultspath, "filterAndTrim_out.rds"))
out <- readRDS(file.path(resultspath, "filterAndTrim_out.rds"))

# Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Save and load error rates models
saveRDS(errF, file.path(resultspath, "errF.rds"))
saveRDS(errR, file.path(resultspath, "errR.rds"))
errF <- readRDS(file.path(resultspath, "errF.rds"))
errR <- readRDS(file.path(resultspath, "errR.rds"))

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

# Save and load inferences
saveRDS(dadaFs, file.path(resultspath, "dadaFs.rds"))
saveRDS(dadaRs, file.path(resultspath, "dadaRs.rds"))
dadaFs <- readRDS(file.path(resultspath, "dadaFs.rds"))
dadaRs <- readRDS(file.path(resultspath, "dadaRs.rds"))

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct Amplicon Sequence Variant (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

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
taxa <- assignTaxonomy(seqtab.nochim, "~/DADA2/tax/silva_nr99_v138.2_toGenus_trainset.fa.gz", minBoot = 50, multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Taxonomy summary
colSums(!is.na(taxa)) # ASVs assigned in each level
round(100 * colSums(!is.na(taxa)) / nrow(taxa),2)

head(sort(table(taxa[, "Phylum"]), decreasing = TRUE), 15) # Main phylums
head(sort(table(taxa[, "Class"]), decreasing=TRUE), 15) # Main classes
head(sort(table(taxa[, "Genus"]), decreasing=TRUE), 20) # Main genus

#######~~~~~~~~~~~~~~~~~~~~~~~~~ Handoff to phyloseq ~~~~~~~~~~~~~~~~~~~~~~~~~#######
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

# Load metadata from sample_data
metadata <- read.csv(
  "~/DADA2.Analysis/metadata/sample_data.csv",
  sep = ";",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)
all(rownames(metadata) == rownames(seqtab.nochim)) # Verify sample names

# Convert variable types
metadata$Replicate <- factor(metadata$Replicate)
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
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(metadata)
)

# Save phyloseq object
saveRDS(ps, file.path(resultspath, "phyloseq_object.rds"))
ps <- readRDS(file.path(resultspath, "phyloseq_object.rds"))

# Rename taxa to short names of ASVs and store the DNA sequences in refseq slot
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Visualize alpha-diversity
plot_richness(
  ps,
  x = "CropStage",
  measures = c("Shannon", "Observed"),
  color = "Year"
)

# Ordinate
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")


plot_ordination(ps.prop, ord.nmds.bray, color="CropStage", title="Bray NMDS") # Plot Bray NMDS

# Relative taxonomic composition
ps.family <- tax_glom(ps, taxrank = "Family")
top20 <- names(sort(taxa_sums(ps.family), decreasing = TRUE))[1:6]
ps.family.rel <- transform_sample_counts(ps.family, function(x) x / sum(x))
ps.family.top20 <- prune_taxa(top20, ps.family.rel)

plot_bar(
  ps.family.top20,
  x = "CropStage",
  fill = "Family"
) +
  facet_wrap(~Year, scales="free_x")
