#####################################################################################
# Project: Coffe Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: load_objects.R
# Autor: Gustavo Goncalves
# Last edit: 29/07/2026
#
# Description:
#    Script to load up objects saved previously on "path".
#
# Current path: ~/results/objects
#
#####################################################################################

# Set up path
path <- "results/objects"; list.files(path)

# Loads
out <- readRDS(file.path(path, "filterAndTrim_out.rds"))
errF <- readRDS(file.path(path, "errF.rds"))
errR <- readRDS(file.path(path, "errR.rds"))
dadaFs <- readRDS(file.path(path, "dadaFs.rds"))
dadaRs <- readRDS(file.path(path, "dadaFs.rds"))
seqtab <- readRDS(file.path(path, "seqtab.rds"))
seqtab.nochim <- readRDS(file.path(path, "seqtab-nochim.rds"))
track <- read.csv(file.path(dirname(path), "read_tracking.csv"),row.names = 1)
taxa <- readRDS(file.path(path, "taxonomy.rds"))
ps <- readRDS(file.path(path, "phyloseq_object.rds"))

# Set ggplot theme
theme_set(theme_bw())
