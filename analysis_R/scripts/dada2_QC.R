#####################################################################################
# Project: Coffee Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: dada2_QC.R
# Author: Gustavo Goncalves
# Last edit: 17/08/2026
#
# Description:
# Quality control report for dada2_asv_inference.R
#
# Required objects:
#   - seqtab.nochim
#   - track
#   - taxa
#
#####################################################################################

# Set up packages
library(dada2); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")

theme_set(theme_bw())


#####################################################################################
# QC calculations
#####################################################################################

# General information
n_samples <- nrow(seqtab.nochim)
n_asvs <- ncol(seqtab.nochim)

reads <- rowSums(seqtab.nochim)

mean_reads <- mean(reads)
median_reads <- median(reads)

# Read retention
track$Retention <- round(
  100 * track$nonchim / track$input,
  2
)

retention_summary <- summary(track$Retention)
mean_retention <- mean(track$Retention)
median_retention <- median(track$Retention)

# ASVs per sample
asvs <- apply(seqtab.nochim > 0, 1, sum)

# ASV abundance
abundance <- colSums(seqtab.nochim)

# ASV lengths
lengths <- nchar(
  getSequences(seqtab.nochim)
)

# ASVs <400 bp
reads_short <- sum(abundance[lengths < 400])
percent_reads_short <- 100 * reads_short / sum(abundance)

# Taxonomy
assigned <- colSums(!is.na(taxa))
assigned_percent <- round(
  100 * assigned / nrow(taxa),
  2
)

phylum <- sort(
  table(taxa[, "Phylum"]),
  decreasing = TRUE
)

family <- sort(
  table(taxa[, "Family"]),
  decreasing = TRUE
)


#####################################################################################
# Visual QC
#####################################################################################

# Read retention
ggplot(
  track,
  aes(y = Retention)
) +
  geom_boxplot()

# ASVs per sample
ggplot(
  data.frame(
    Sample = names(asvs),
    ASVs = asvs
  ),
  aes(Sample, ASVs)
) +
  geom_col() +
  coord_flip()

# ASV abundance
hist(
  log10(abundance),
  breaks = 50
)

# ASV lengths
hist(
  lengths,
  breaks = 200
)

# Most abundant phyla
phylum_top20 <- head(phylum, 20)

ggplot(
  data.frame(
    Phylum = names(phylum_top20),
    Count = as.numeric(phylum_top20)
  ),
  aes(x = reorder(Phylum, Count), y = Count)
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Phylum",
    y = "Number of ASVs"
  )

# Most abundant families
family_top20 <- head(family, 20)

ggplot(
  data.frame(
    Family = names(family_top20),
    Count = as.numeric(family_top20)
  ),
  aes(x = reorder(Family, Count), y = Count)
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Family",
    y = "Number of ASVs"
  )


#####################################################################################
# Detect possible problematic samples
#####################################################################################

cat("Samples with <30% retention:\n")

track[
  track$Retention < 30,
  c("input", "nonchim", "Retention")
]

cat("\nSamples with <1000 ASVs:\n")

asvs[asvs < 1000]

cat("\nSamples with <10000 reads:\n")

reads[reads < 10000]


#####################################################################################
# Save QC summary report
#####################################################################################

sink("results/dada2-QC/dada2_QC_summary.txt")

cat("=====================================\n")
cat("DADA2 QC SUMMARY REPORT\n")
cat("Coffee Rhizosphere Microbiome - 16S rRNA\n")
cat("=====================================\n\n")

cat("Generated:", format(Sys.time()), "\n\n")


# General information
cat("-------------------------------------\n")
cat("GENERAL INFORMATION\n")
cat("-------------------------------------\n\n")

cat("Number of samples:", n_samples, "\n")
cat("Total ASVs:", n_asvs, "\n")
cat("Total reads:", sum(reads), "\n")
cat("Mean reads/sample:", round(mean_reads), "\n")
cat("Median reads/sample:", round(median_reads), "\n\n")


# Read retention
cat("-------------------------------------\n")
cat("READ RETENTION\n")
cat("-------------------------------------\n\n")

cat("Retention summary (%):\n")
print(retention_summary)

cat("\nMean retention (%):",
    round(mean_retention, 2), "\n")

cat("Median retention (%):",
    round(median_retention, 2), "\n\n")


# Problematic samples
cat("-------------------------------------\n")
cat("POSSIBLE PROBLEMATIC SAMPLES\n")
cat("-------------------------------------\n\n")

cat("Samples with <30% retention:\n")
print(
  track[
    track$Retention < 30,
    c("input", "nonchim", "Retention")
  ]
)

cat("\nSamples with <1000 ASVs:\n")
print(asvs[asvs < 1000])

cat("\nSamples with <10000 reads:\n")
print(reads[reads < 10000])


# ASV summary
cat("\n-------------------------------------\n")
cat("ASV SUMMARY\n")
cat("-------------------------------------\n\n")

cat("ASVs per sample:\n")
print(
  data.frame(
    Sample = names(asvs),
    ASVs = asvs
  )
)

cat("\nASV distribution:\n")
print(summary(asvs))

cat("\nASV abundance distribution:\n")
print(summary(abundance))

cat("\nSingleton ASVs:",
    sum(abundance == 1), "\n")

cat(
  "ASVs represented by singletons (%):",
  round(
    100 * sum(abundance == 1) / length(abundance),
    2
  ),
  "\n\n"
)


# ASV lengths
cat("-------------------------------------\n")
cat("ASV LENGTHS\n")
cat("-------------------------------------\n\n")

cat("Length distribution:\n")
print(table(lengths))

cat("\nReads from ASVs <400 bp:",
    reads_short, "\n")

cat(
  "Percentage of reads from ASVs <400 bp:",
  round(percent_reads_short, 2),
  "%\n"
)


# Taxonomy
cat("\n-------------------------------------\n")
cat("TAXONOMIC ASSIGNMENT\n")
cat("-------------------------------------\n\n")

cat("Assigned taxa per rank:\n")
print(assigned)

cat("\nPercentage assigned:\n")
print(assigned_percent)

cat("\nMost abundant Phyla:\n")
print(head(phylum, 20))

cat("\nMost abundant Families:\n")
print(head(family, 20))

sink()

#####################################################################################
# End QC
#####################################################################################
