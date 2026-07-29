#####################################################################################
# Project: Coffe Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: dada2_QC.R
# Autor: Gustavo Goncalves
# Last edit: 17/07/2026
#
# Description:
# Quality control report for dada2_asv_inference.R:
#
#####################################################################################

# Set up packages
library(dada2); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")

theme_set(theme_bw())

# Input summary
cat("=====================================\n")
cat("DADA2 QC REPORT\n")
cat("=====================================\n\n")
cat("Samples processed:", nrow(seqtab.nochim), "\n")
cat("Total ASVs:", ncol(seqtab.nochim), "\n\n")
cat("Mean reads/sample:", round(mean(rowSums(seqtab.nochim))), "\n")
cat("Median reads/sample:", round(median(rowSums(seqtab.nochim))), "\n")

# Tracking read retention
track$Retention <- round(
  100 * track$nonchim / track$input,
  2
)
summary(track$Retention)
cat("Mean retention (%):", round(mean(track$Retention),2), "\n")

cat("Median retention (%):", round(median(track$Retention),2), "\n")
ggplot(track,
       aes(y = Retention)) +
  geom_boxplot()

# ASV per sample
asvs <- apply(seqtab.nochim > 0, 1, sum)
summary(asvs)

ggplot(data.frame(Sample = names(asvs), ASVs = asvs), aes(Sample, ASVs)
) +
  geom_col() +
  coord_flip()

# ASVs abundance
abundance <- colSums(seqtab.nochim)

summary(abundance)
hist(
  log10(abundance),
  breaks = 50
)

# ASVs lengths
lengths <- nchar(
  getSequences(seqtab.nochim)
)

table(lengths)

hist(
  lengths,
  breaks = 200
)
sum(abundance[lengths < 400])
sum(abundance)

100 * sum(abundance[lengths < 400]) / sum(abundance)

 # Taxonomy
assigned <- colSums(!is.na(taxa))
assigned
round(100 * assigned / nrow(taxa), 2)

phylum <- sort(
  table(
    taxa[, "Phylum"]),
  decreasing = TRUE)

head(phylum,20)
barplot(
  head(phylum,20),
  las=2
)

family <- sort(
  table(
    taxa[, "Family"]),
  decreasing = TRUE)

head(family,20)
barplot(
  head(family,20),
  las=2
)

# Phyloseq object
ps

# Detect possible problematic samples

cat("Samples with <30% retention:\n")
track[track$Retention < 30,
      c("input","nonchim","Retention")]

cat("Samples with <1000 ASVs:\n")
asvs[asvs < 1000]

cat("Samples with <10000 reads:\n")
reads <- rowSums(seqtab.nochim)
reads[reads < 10000]

#####################################################################################
# Save QC summary report
#####################################################################################

sink("dada2_QC_summary.txt")

cat("=====================================\n")
cat("DADA2 QC SUMMARY REPORT\n")
cat("Coffee Rhizosphere Microbiome - 16S rRNA\n")
cat("=====================================\n\n")

cat("Generated:", format(Sys.time()), "\n\n")


# General information
cat("-------------------------------------\n")
cat("GENERAL INFORMATION\n")
cat("-------------------------------------\n\n")

cat("Number of samples:", nrow(seqtab.nochim), "\n")
cat("Total ASVs:", ncol(seqtab.nochim), "\n")

cat("Total reads:", sum(seqtab.nochim), "\n")
cat("Mean reads/sample:",
    round(mean(rowSums(seqtab.nochim))), "\n")
cat("Median reads/sample:",
    round(median(rowSums(seqtab.nochim))), "\n\n")


# Read retention
cat("-------------------------------------\n")
cat("READ RETENTION\n")
cat("-------------------------------------\n\n")

print(track)

cat("\nRetention summary (%):\n")
print(summary(track$Retention))

cat("\nMean retention (%):",
    round(mean(track$Retention),2), "\n")

cat("Median retention (%):",
    round(median(track$Retention),2), "\n\n")


# Low quality samples
cat("-------------------------------------\n")
cat("POSSIBLE PROBLEMATIC SAMPLES\n")
cat("-------------------------------------\n\n")

cat("Samples with <30% retention:\n")
print(track[track$Retention < 30,
            c("input","nonchim","Retention")])

cat("\nSamples with <1000 ASVs:\n")
print(asvs[asvs < 1000])

cat("\nSamples with <10000 reads:\n")
print(reads[reads < 10000])


# ASV summary
cat("\n\n-------------------------------------\n")
cat("ASV SUMMARY\n")
cat("-------------------------------------\n\n")

cat("ASVs per sample:\n")
print(data.frame(
  Sample = names(asvs),
  ASVs = asvs
))

cat("\nASV distribution:\n")
print(summary(asvs))


cat("\nASV abundance distribution:\n")
print(summary(abundance))


cat("\nSingleton ASVs:",
    sum(abundance == 1), "\n")

cat("ASVs represented by singletons (%):",
    round(100 * sum(abundance == 1) /
            length(abundance),2),
    "\n\n")


# ASV length
cat("-------------------------------------\n")
cat("ASV LENGTHS\n")
cat("-------------------------------------\n\n")

print(table(lengths))

cat("\nReads from ASVs <400 bp:\n")
print(sum(abundance[lengths < 400]))

cat("Percentage of reads from ASVs <400 bp:\n")
print(round(
  100 * sum(abundance[lengths < 400]) /
    sum(abundance),
  2
))


# Taxonomy assignment
cat("\n\n-------------------------------------\n")
cat("TAXONOMIC ASSIGNMENT\n")
cat("-------------------------------------\n\n")

assigned <- colSums(!is.na(taxa))

cat("Assigned taxa per rank:\n")
print(assigned)

cat("\nPercentage assigned:\n")
print(round(
  100 * assigned / nrow(taxa),
  2
))


cat("\nMost abundant Phyla:\n")
print(head(phylum,20))


cat("\nMost abundant Families:\n")
print(head(family,20))


# Phyloseq summary
cat("\n\n-------------------------------------\n")
cat("PHYLOSEQ OBJECT\n")
cat("-------------------------------------\n\n")

print(ps)


cat("\nNumber of taxa:")
print(ntaxa(ps))

cat("\nNumber of samples:")
print(nsamples(ps))

sink()

#####################################################################################
# End QC summary
#####################################################################################
