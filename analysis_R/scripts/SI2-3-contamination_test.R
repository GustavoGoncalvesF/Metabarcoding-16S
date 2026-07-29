#####################################################################################
# Script to test SI2 and SI3 contamination
#####################################################################################

alpha <- estimate_richness(ps)
alpha[c("SI2","SI3"),]

ps.bac <- subset_taxa(ps, Family == "Bacillaceae")
otu_table(ps.bac)

bac <- as(otu_table(ps.bac), "matrix")
bac_SI2 <- sort(bac["SI2", ], decreasing = TRUE)
head(bac_SI2, 20)
bac_SI3 <- sort(bac["SI3", ], decreasing = TRUE)
head(bac_SI3, 20)

tax_table(ps)[c("ASV1","ASV4","ASV37","ASV81","ASV162"), ]
c("ASV1","ASV4","ASV37","ASV81","ASV162")
otu <- as(otu_table(ps), "matrix")
otu.sub <- otu[, c("ASV1","ASV4","ASV37","ASV81","ASV162")]
otu.rel <- otu.sub / rowSums(otu)

library(reshape2)

df <- melt(otu.rel)

colnames(df) <- c("Sample","ASV","RelativeAbundance")

ggplot(df,
       aes(ASV,
           Sample,
           fill = RelativeAbundance)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw()
