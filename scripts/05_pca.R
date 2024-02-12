# Load the required packages ----------------------------------------------

# For determining the project root
library(here)
# For transforming, representing, and plotting data
library(tidyverse)
library(ggrepel)
library(plotly)
# For data analysis
library(vcfR)
library(adegenet)


# Set scientific notation options and default plotting theme --------------

options(scipen = 1)
options(pillar.sigfig = 4)

ggplot2::theme_set(ggplot2::theme_classic())


# File paths --------------------------------------------------------------

dta_accessions_file_path <- here::here(
  "data",
  "accessions.tsv"
)

dta_vcf_file_path <- here::here(
  "data",
  "final_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz"
)


# Read in the data --------------------------------------------------------

# Read in the accessions data
dta_accessions <- readr::read_tsv(
  file = dta_accessions_file_path,
  col_types = "ffnn"
)

# Read in the VCF
dta_vcf <- vcfR::read.vcfR(file = dta_vcf_file_path)


# Principle Component Analysis (PCA) --------------------------------------

# Convert VCF to genind format, scale for PCA and perform PCA
pca <- vcfR::vcfR2genind(dta_vcf) |>
  adegenet::scaleGen(NA.method = "mean") |>
  ade4::dudi.pca(scale = FALSE, scannf = FALSE, nf = 10)

pca_axis_df <- pca$li
# Add a new column "ind" to pca_axis using the accession names
pca_axis_df$ind <- dta_accessions$sample
# Add the population information to pca_axis_df
pca_axis_df$Population <- dta_accessions$pop
# Convert to a tibble
pca_axis_df <- tibble::as_tibble(pca_axis_df)

# Percentage variance explained
pca_eigenval <- pca$eig[1:10]
pve_df <- tibble::as_tibble(data.frame(
  PC = 1:10,
  pve = pca_eigenval / sum(pca_eigenval) * 100
))
# Calculate the cumulative sum of the percentage variance explained
cumsum(pve_df$pve)

# Plot the percentage of variance explained
ggplot2::ggplot(pve_df, ggplot2::aes(PC, pve)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::labs(
    title = "PCA - A. thaliana accessions from Spain (ESP) and Sweden (SWE)",
    y = "Percentage of variance explained"
  )

# Plot PC1 and PC2
ggplot2::ggplot(
  pca_axis_df,
  ggplot2::aes(Axis1, Axis2, col = Population, label = ind)
) +
  ggplot2::geom_point(size = 2) +
  ggrepel::geom_text_repel(
    show.legend = FALSE,
    size = 2.5,
    fontface = "bold",
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  ggplot2::scale_colour_manual(values = c("blue", "red")) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    x = paste0("PC1 (", signif(pve_df$pve[1], 3), "%)"),
    y = paste0("PC2 (", signif(pve_df$pve[2], 3), "%)"),
    title = "PCA - A. thaliana accessions from Spain (ESP)\nand Sweden (SWE)",
  ) +
  ggplot2::theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

# PCA visualization using interactive 3D plot
plotly::plot_ly(
  pca_axis_df,
  x = ~Axis1,
  y = ~Axis2,
  z = ~Axis3,
  color = dta_accessions$pop,
  colors = c("blue", "red")
) |>
  plotly::add_markers(size = 12) |>
  plotly::layout(
    title = "PCA - A. thaliana accessions from Spain (ESP) and Sweden (SWE)",
    scene = list(bgcolor = "#e5ecf6"),
    plot_bgcolor = "rgb(0, 0, 0,0)",
    paper_bgcolor = "rgb(0, 0, 0,0)",
    margin = list(l = 80, r = 80, b = 80, t = 80),
    legend = list(bgcolor = "rgba(0,0,0,0)")
  )
