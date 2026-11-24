#####################################################################################
# Project: 16S rRNA Metabarcoding
# Script: dada2_read_tracking.R
# Autor: Gustavo Goncalves
# Last edit: 02/11/2026
#
# Description:
#   Auxiliary script to produce the appropriate read tracking originating from
#   01_dada2_asv_inference.R;
#
#####################################################################################

# Set up packages
library(openxlsx); packageVersion("openxlsx")

generate_read_tracking <- function(
    out, dadaFs, dadaRs, mergers,
    seqtab.nochim, sample.names, results_path
) {

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required.")
  }

  getN <- function(x) sum(dada2::getUniques(x))

  # Read counts
  track <- data.frame(
    SampleID = sample.names,
    input = as.numeric(out[, 1]),
    filtered = as.numeric(out[, 2]),
    denoisedF = as.numeric(sapply(dadaFs, getN)),
    denoisedR = as.numeric(sapply(dadaRs, getN)),
    merged = as.numeric(sapply(mergers, getN)),
    nonchim = as.numeric(rowSums(seqtab.nochim)),
    row.names = NULL
  )

  # Retention percentages
  pct <- function(a, b) {
    ifelse(b > 0, 100 * a / b, NA_real_)
  }

  retention <- data.frame(
    SampleID = track$SampleID,
    Filtering = pct(track$filtered, track$input),
    DenoisedF = pct(track$denoisedF, track$filtered),
    DenoisedR = pct(track$denoisedR, track$filtered),
    Merging = pct(track$merged, track$filtered),
    Nonchim = pct(track$nonchim, track$merged),
    Total = pct(track$nonchim, track$input)
  )

  # Absolute losses
  losses <- data.frame(
    SampleID = track$SampleID,
    Filtering = track$input - track$filtered,
    DenoisedF = track$filtered - track$denoisedF,
    DenoisedR = track$filtered - track$denoisedR,
    Merging = track$filtered - track$merged,
    Chimeras = track$merged - track$nonchim,
    Total = track$input - track$nonchim
  )

  # Quality flags
  track$QC_Flag <- ifelse(
    track$nonchim < 10000, "LOW_DEPTH",
    ifelse(retention$Total < 30, "LOW_RETENTION", "OK")
  )

  # Complete CSV
  complete <- cbind(track, retention[-1], losses[-1])

  dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    complete,
    file.path(results_path, "read_tracking.csv"),
    row.names = FALSE
  )

  # Excel workbook
  wb <- openxlsx::createWorkbook()

  sheets <- list(
    "01_Read_Counts" = track,
    "02_Retention" = retention,
    "03_Read_Losses" = losses
  )

  for (sheet in names(sheets)) {
    openxlsx::addWorksheet(wb, sheet)

    openxlsx::writeDataTable(
      wb, sheet, sheets[[sheet]],
      tableStyle = "TableStyleMedium2"
    )

    openxlsx::freezePane(
      wb, sheet, firstActiveRow = 2
    )

    openxlsx::setColWidths(
      wb, sheet,
      cols = seq_len(ncol(sheets[[sheet]])),
      widths = "auto"
    )
  }

  # Highlight low final read counts
  openxlsx::conditionalFormatting(
    wb, "01_Read_Counts",
    cols = 7,
    rows = 2:(nrow(track) + 1),
    rule = "<10000",
    style = openxlsx::createStyle(
      bgFill = "#F4CCCC"
    )
  )

  # Highlight low total retention
  openxlsx::conditionalFormatting(
    wb, "02_Retention",
    cols = ncol(retention),
    rows = 2:(nrow(track) + 1),
    rule = "<30",
    style = openxlsx::createStyle(
      bgFill = "#F4CCCC"
    )
  )

  openxlsx::saveWorkbook(
    wb,
    file.path(results_path, "read_tracking.xlsx"),
    overwrite = TRUE
  )

  message("Read tracking generated successfully.")

  invisible(list(
    counts = track,
    retention = retention,
    losses = losses
  ))
}
