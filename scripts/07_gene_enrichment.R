# Load the required packages ----------------------------------------------

# For determining the project root
library(here)
# For transforming, representing, and plotting data
library(tidyverse)
# For data analysis
library(topGO)


# Set scientific notation options and default plotting theme --------------

options(scipen = 1)
options(pillar.sigfig = 4)

ggplot2::theme_set(ggplot2::theme_classic())


# File paths --------------------------------------------------------------

dta_gene_list_file_path <- here::here(
  "data",
  "GO_input.txt"
)

dta_enriched_go_mf_file_path <- here::here(
  "data",
  "enriched_go_mf.csv"
)

dta_enriched_go_bp_file_path <- here::here(
  "data",
  "enriched_go_bp.csv"
)

dta_enriched_go_cc_file_path <- here::here(
  "data",
  "enriched_go_cc.csv"
)


# Read in the data --------------------------------------------------------

# Read the gene expression and significance values
gene_list <- read.table(dta_gene_list_file_path, header = TRUE)


# Gene enrichment analysis ------------------------------------------------

# Create a universe (list of all genes)
# for enrichments based on significance

# Select the third column with p-values
univ <- gene_list[, 3]
# Give names to genes from the first column
names(univ) <- gene_list[, 1]
# Ignore NA values
univ <- univ[!is.na(univ)]


# MF enrichment -----------------------------------------------------------

# Prepare topGO data
tGOdataMF <- new(
  "topGOdata",
  description = "Simple session",
  ontology = "MF",
  geneSel = function(score) score < 0.05,
  allGenes = univ,
  nodeSize = 3,
  mapping = "org.At.tair.db",
  annot = annFUN.org
)

# Run Gene Enrichment Analysis
results_mf_ks <- topGO::runTest(tGOdataMF, algorithm = "elim", statistic = "ks")

# Generate a table of enriched GO terms
goEnrichmentMF <- topGO::GenTable(
  tGOdataMF,
  KS = results_mf_ks,
  orderBy = "KS",
  topNodes = 50
)
goEnrichmentMF$KS <- as.numeric(goEnrichmentMF$KS)
goEnrichmentMF <- goEnrichmentMF[goEnrichmentMF$KS < 0.05, ]
goEnrichmentMF <- goEnrichmentMF[, c("GO.ID", "Term", "KS")]
goEnrichmentMF$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichmentMF$Term)
goEnrichmentMF$Term <- gsub("\\.\\.\\.$", "", goEnrichmentMF$Term)
goEnrichmentMF$Term <- paste(
  goEnrichmentMF$GO.ID,
  goEnrichmentMF$Term,
  sep = ", ")

# Remove NAs if any exist
if (any(is.na(goEnrichmentMF$KS))) {
  goEnrichmentMF <- goEnrichmentMF[!is.na(goEnrichmentMF$KS), ]
}

# Save the enriched GO terms table to a CSV file
write.csv(goEnrichmentMF, dta_enriched_go_mf_file_path)

# Select the top N terms
ntop <- 10
ggdata_mf <- goEnrichmentMF[1:ntop, ]

# Fix order
ggdata_mf$Term <- factor(ggdata_mf$Term, levels = rev(ggdata_mf$Term))

# Create a plot
ggplot2::ggplot(
  ggdata_mf,
  ggplot2::aes(x = Term, y = -log10(KS), size = -log10(KS), fill = -log10(KS))
) +
  ggplot2::expand_limits(y = 1) +
  ggplot2::geom_point(shape = 21) +
  ggplot2::scale_size(range = c(2.5, 12.5)) +
  ggplot2::scale_fill_continuous(low = "royalblue", high = "red4") +
  ggplot2::labs(
    title = "Pathways Enriched in Genes with High Fst Values - MF",
    subtitle = "Top 10 Terms Ordered by Kolmogorov-Smirnov p-Value",
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001",
    x = NULL,
    y = "Enrichment score"
  ) +
  ggplot2::geom_hline(
    yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    linetype = c("dotted", "longdash", "solid"),
    colour = c("black", "black", "black"),
    linewidth = c(0.5, 1.5, 3)
  ) +
  ggplot2::theme(
    legend.position = "right",
    legend.background = ggplot2::element_rect(),
    plot.title = ggplot2::element_text(size = 18, face = "bold", vjust = 1),
    plot.subtitle = ggplot2::element_text(size = 16, face = "bold", vjust = 1),
    plot.caption = ggplot2::element_text(size = 12, vjust = 1),
    axis.text.x = ggplot2::element_text(size = 12, hjust = 1.10),
    axis.text.y = ggplot2::element_text(
      size = 12,
      face = "bold",
      vjust = 0.5),
    axis.title = ggplot2::element_text(size = 12),
    axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
    axis.line = ggplot2::element_line(colour = "black"),
    legend.key = ggplot2::element_blank(),
    legend.key.size = ggplot2::unit(1, "cm"),
    legend.text = ggplot2::element_text(size = 16),
    title = ggplot2::element_text(size = 12)
  ) +
  ggplot2::coord_flip()


# BP enrichment -----------------------------------------------------------

# Prepare topGO data
tGOdataBP <- new(
  "topGOdata",
  description = "Simple session",
  ontology = "BP",
  geneSel = function(score) score < 0.05,
  allGenes = univ,
  nodeSize = 3,
  mapping = "org.At.tair.db",
  annot = annFUN.org
)

# Run Gene Enrichment Analysis
results_bp_ks <- topGO::runTest(tGOdataBP, algorithm = "elim", statistic = "ks")

# Generate a table of enriched GO terms
goEnrichmentBP <- topGO::GenTable(
  tGOdataBP,
  KS = results_bp_ks,
  orderBy = "KS",
  topNodes = 50
)
goEnrichmentBP$KS <- as.numeric(goEnrichmentBP$KS)
goEnrichmentBP <- goEnrichmentBP[goEnrichmentBP$KS < 0.05, ]
goEnrichmentBP <- goEnrichmentBP[, c("GO.ID", "Term", "KS")]
goEnrichmentBP$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichmentBP$Term)
goEnrichmentBP$Term <- gsub("\\.\\.\\.$", "", goEnrichmentBP$Term)
goEnrichmentBP$Term <- paste(
  goEnrichmentBP$GO.ID,
  goEnrichmentBP$Term,
  sep = ", ")

# Remove NAs if any exist
if (any(is.na(goEnrichmentBP$KS))) {
  goEnrichmentBP <- goEnrichmentBP[!is.na(goEnrichmentBP$KS), ]
}

# Save the enriched GO terms table to a CSV file
write.csv(goEnrichmentBP, dta_enriched_go_bp_file_path)

# Select the top N terms
ntop <- 10
ggdata_bp <- goEnrichmentBP[1:ntop, ]

# Fix order
ggdata_bp$Term <- factor(ggdata_bp$Term, levels = rev(ggdata_bp$Term))

# Create a plot
ggplot2::ggplot(
  ggdata_bp,
  ggplot2::aes(x = Term, y = -log10(KS), size = -log10(KS), fill = -log10(KS))
) +
  ggplot2::expand_limits(y = 1) +
  ggplot2::geom_point(shape = 21) +
  ggplot2::scale_size(range = c(2.5, 12.5)) +
  ggplot2::scale_fill_continuous(low = "royalblue", high = "red4") +
  ggplot2::labs(
    title = "Pathways Enriched in Genes with High Fst Values - BP",
    subtitle = "Top 10 Terms Ordered by Kolmogorov-Smirnov p-Value",
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001",
    x = NULL,
    y = "Enrichment score"
  ) +
  ggplot2::geom_hline(
    yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    linetype = c("dotted", "longdash", "solid"),
    colour = c("black", "black", "black"),
    linewidth = c(0.5, 1.5, 3)
  ) +
  ggplot2::theme(
    legend.position = "right",
    legend.background = ggplot2::element_rect(),
    plot.title = ggplot2::element_text(size = 18, face = "bold", vjust = 1),
    plot.subtitle = ggplot2::element_text(size = 16, face = "bold", vjust = 1),
    plot.caption = ggplot2::element_text(size = 12, vjust = 1),
    axis.text.x = ggplot2::element_text(size = 12, hjust = 1.10),
    axis.text.y = ggplot2::element_text(
      size = 12,
      face = "bold",
      vjust = 0.5),
    axis.title = ggplot2::element_text(size = 12),
    axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
    axis.line = ggplot2::element_line(colour = "black"),
    legend.key = ggplot2::element_blank(),
    legend.key.size = ggplot2::unit(1, "cm"),
    legend.text = ggplot2::element_text(size = 16),
    title = ggplot2::element_text(size = 12)
  ) +
  ggplot2::coord_flip()


# CC enrichment -----------------------------------------------------------

# Prepare topGO data
tGOdataCC <- new(
  "topGOdata",
  description = "Simple session",
  ontology = "CC",
  geneSel = function(score) score < 0.05,
  allGenes = univ,
  nodeSize = 3,
  mapping = "org.At.tair.db",
  annot = annFUN.org
)

# Run Gene Enrichment Analysis
results_cc_ks <- topGO::runTest(tGOdataCC, algorithm = "elim", statistic = "ks")

# Generate a table of enriched GO terms
goEnrichmentCC <- topGO::GenTable(
  tGOdataCC,
  KS = results_cc_ks,
  orderBy = "KS",
  topNodes = 50
)
goEnrichmentCC$KS <- as.numeric(goEnrichmentCC$KS)
goEnrichmentCC <- goEnrichmentCC[goEnrichmentCC$KS < 0.05, ]
goEnrichmentCC <- goEnrichmentCC[, c("GO.ID", "Term", "KS")]
goEnrichmentCC$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichmentCC$Term)
goEnrichmentCC$Term <- gsub("\\.\\.\\.$", "", goEnrichmentCC$Term)
goEnrichmentCC$Term <- paste(
  goEnrichmentCC$GO.ID,
  goEnrichmentCC$Term,
  sep = ", ")

# Remove NAs if any exist
if (any(is.na(goEnrichmentCC$KS))) {
  goEnrichmentCC <- goEnrichmentCC[!is.na(goEnrichmentCC$KS), ]
}

# Save the enriched GO terms table to a CSV file
write.csv(goEnrichmentCC, dta_enriched_go_cc_file_path)

# Select the top N terms
ntop <- 10
ggdata_cc <- goEnrichmentCC[1:ntop, ]

# Fix order
ggdata_cc$Term <- factor(ggdata_cc$Term, levels = rev(ggdata_cc$Term))

# Create a plot
ggplot2::ggplot(
  ggdata_cc,
  ggplot2::aes(x = Term, y = -log10(KS), size = -log10(KS), fill = -log10(KS))
) +
  ggplot2::expand_limits(y = 1) +
  ggplot2::geom_point(shape = 21) +
  ggplot2::scale_size(range = c(2.5, 12.5)) +
  ggplot2::scale_fill_continuous(low = "royalblue", high = "red4") +
  ggplot2::labs(
    title = "Pathways Enriched in Genes with High Fst Values - CC",
    subtitle = "Top 10 Terms Ordered by Kolmogorov-Smirnov p-Value",
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001",
    x = NULL,
    y = "Enrichment score"
  ) +
  ggplot2::geom_hline(
    yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    linetype = c("dotted", "longdash", "solid"),
    colour = c("black", "black", "black"),
    linewidth = c(0.5, 1.5, 3)
  ) +
  ggplot2::theme(
    legend.position = "right",
    legend.background = ggplot2::element_rect(),
    plot.title = ggplot2::element_text(size = 18, face = "bold", vjust = 1),
    plot.subtitle = ggplot2::element_text(size = 16, face = "bold", vjust = 1),
    plot.caption = ggplot2::element_text(size = 12, vjust = 1),
    axis.text.x = ggplot2::element_text(size = 12, hjust = 1.10),
    axis.text.y = ggplot2::element_text(
      size = 12,
      face = "bold",
      vjust = 0.5),
    axis.title = ggplot2::element_text(size = 12),
    axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
    axis.line = ggplot2::element_line(colour = "black"),
    legend.key = ggplot2::element_blank(),
    legend.key.size = ggplot2::unit(1, "cm"),
    legend.text = ggplot2::element_text(size = 16),
    title = ggplot2::element_text(size = 12)
  ) +
  ggplot2::coord_flip()
