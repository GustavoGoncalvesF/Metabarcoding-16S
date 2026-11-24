#####################################################################################
# Pairwise Bray-Curtis dissimilarity
#####################################################################################

# Extract Bray-Curtis distance matrix
bray.mat <- as.matrix(bray.balanced)

# Extract sample pairs
pairs <- t(combn(rownames(bray.mat), 2))

# Create pairwise comparison table
pairwise <- data.frame(
  Sample1 = pairs[, 1],
  Sample2 = pairs[, 2]
)

# Add Bray-Curtis distances
pairwise$Bray <- mapply(
  function(a, b) bray.mat[a, b],
  pairwise$Sample1,
  pairwise$Sample2
)

# Identify shared species
pairwise$SameSpecies <-
  metadata.balanced[pairwise$Sample1, "Species"] ==
  metadata.balanced[pairwise$Sample2, "Species"]

# Identify shared origin
pairwise$SameOrigin <-
  metadata.balanced[pairwise$Sample1, "Origin"] ==
  metadata.balanced[pairwise$Sample2, "Origin"]

# Classify comparisons
pairwise$Comparison <- with(
  pairwise,
  ifelse(
    SameSpecies & SameOrigin,
    "Same species / Same origin",
    ifelse(
      SameSpecies & !SameOrigin,
      "Same species / Different origin",
      ifelse(
        !SameSpecies & SameOrigin,
        "Different species / Same origin",
        "Different species / Different origin"
      )
    )
  )
)

# Set comparison order
pairwise$Comparison <- factor(
  pairwise$Comparison,
  levels = c(
    "Same species / Same origin",
    "Same species / Different origin",
    "Different species / Same origin",
    "Different species / Different origin"
  )
)

# Summary statistics
aggregate(
  Bray ~ Comparison,
  data = pairwise,
  FUN = function(x) c(
    Mean = mean(x),
    Median = median(x),
    SD = sd(x)
  )
)


#####################################################################################
# Visualization
#####################################################################################

p.pairwise <- ggplot(
  pairwise,
  aes(
    x = Comparison,
    y = Bray,
    fill = Comparison
  )
) +
  geom_boxplot(alpha = 0.6) +
  labs(
    x = "Community comparison",
    y = "Bray-Curtis dissimilarity"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    )
  ) +
  guides(fill = "none")

print(p.pairwise)

ggsave(
  file.path(beta_path, "10_pairwise_bray_curtis.png"),
  plot = p.pairwise,
  width = 11,
  height = 6,
  dpi = 300
)


#####################################################################################
# Refined pairwise Bray-Curtis comparisons
#####################################################################################

# Extract Bray-Curtis distance matrix
bray.mat <- as.matrix(bray.balanced)

# Generate sample pairs
pairs <- t(combn(rownames(bray.mat), 2))

pairwise.refined <- data.frame(
  Sample1 = pairs[, 1],
  Sample2 = pairs[, 2]
)

# Extract Bray-Curtis dissimilarities
pairwise.refined$Bray <- mapply(
  function(a, b) bray.mat[a, b],
  pairwise.refined$Sample1,
  pairwise.refined$Sample2
)

# Identify shared accession
pairwise.refined$SameAccession <-
  metadata.balanced[pairwise.refined$Sample1, "Accession"] ==
  metadata.balanced[pairwise.refined$Sample2, "Accession"]

# Identify shared species
pairwise.refined$SameSpecies <-
  metadata.balanced[pairwise.refined$Sample1, "Species"] ==
  metadata.balanced[pairwise.refined$Sample2, "Species"]

# Identify shared origin
pairwise.refined$SameOrigin <-
  metadata.balanced[pairwise.refined$Sample1, "Origin"] ==
  metadata.balanced[pairwise.refined$Sample2, "Origin"]


#####################################################################################
# Classify comparisons
#####################################################################################

pairwise.refined$Comparison <- with(
  pairwise.refined,
  ifelse(
    SameAccession & SameOrigin,
    "Same accession / Same origin",
    ifelse(
      SameAccession & !SameOrigin,
      "Same accession / Different origin",
      ifelse(
        SameSpecies & SameOrigin,
        "Different accession / Same species / Same origin",
        ifelse(
          SameSpecies & !SameOrigin,
          "Different accession / Same species / Different origin",
          ifelse(
            SameOrigin,
            "Different species / Same origin",
            "Different species / Different origin"
          )
        )
      )
    )
  )
)

# Set comparison order
pairwise.refined$Comparison <- factor(
  pairwise.refined$Comparison,
  levels = c(
    "Same accession / Same origin",
    "Same accession / Different origin",
    "Different accession / Same species / Same origin",
    "Different accession / Same species / Different origin",
    "Different species / Same origin",
    "Different species / Different origin"
  )
)


#####################################################################################
# Summary statistics
#####################################################################################

summary.refined <- aggregate(
  Bray ~ Comparison,
  data = pairwise.refined,
  FUN = function(x) c(
    N = length(x),
    Mean = mean(x),
    Median = median(x),
    SD = sd(x)
  )
)

print(summary.refined)


#####################################################################################
# Visualization
#####################################################################################

p.refined <- ggplot(
  pairwise.refined,
  aes(
    x = Comparison,
    y = Bray,
    fill = Comparison
  )
) +
  geom_boxplot(alpha = 0.6) +
  labs(
    x = "Community comparison",
    y = "Bray-Curtis dissimilarity"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  ) +
  guides(fill = "none")

print(p.refined)

ggsave(
  file.path(beta_path, "11_pairwise_bray_refined.png"),
  plot = p.refined,
  width = 13,
  height = 7,
  dpi = 300
)

