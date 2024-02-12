# Load the required packages ----------------------------------------------

# For determining the project root
library(here)
# For transforming, representing, and plotting data
library(tidyverse)
# For data analysis
library(VariantAnnotation)
library(PopGenome)


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


# Creates a tibble with the start and end positions of the chromosomes ----

get_chrom_ranges <- function(vcf) {
  chrom_ranges <- tibble::tibble(
    Chromosome = character(),
    `From position` = numeric(),
    `To position` = numeric()
  )
  for (chrom in GenomeInfoDb::seqlevels(vcf)) {
    chrom_pos <- BiocGenerics::start(vcf[GenomeInfoDb::seqnames(vcf) == chrom])
    chrom_ranges <- chrom_ranges |>
      tibble::add_row(
        Chromosome = chrom,
        `From position` = min(chrom_pos),
        `To position` = max(chrom_pos)
      )
  }
  chrom_ranges
}


# Create a tibble with the start and end positions of the chromosones -----

chrom_ranges <- get_chrom_ranges(
  VariantAnnotation::readVcf(file = dta_vcf_file_path)
)


# Extract the population names --------------------------------------------

pop_names <- levels(dta_accessions$pop)


# Create a list with entries for the accessions from Spain and Sweden -----

populations <- split(dta_accessions$sample, dta_accessions$pop)


# Analyze all chromosomes (separately) ------------------------------------

for (chrom in seq_along(chrom_ranges$Chromosome)) {
  from_pos <- dplyr::pull(chrom_ranges[chrom, as.symbol("From position")])
  to_pos <- dplyr::pull(chrom_ranges[chrom, as.symbol("To position")])

  # Read in the VCF file using PopGenome
  dta_chr <- PopGenome::readVCF(
    filename = dta_vcf_file_path,
    numcols = 10000,
    tid = as.character(chrom),
    frompos = from_pos,
    topos = to_pos,
    include.unknown = TRUE
  )

  # Set the populations data
  dta_chr <- PopGenome::set.populations(dta_chr, populations, diploid = TRUE)

  # Print summary of the data
  print(PopGenome::get.sum.data(dta_chr))

  # Create a sliding window dataset
  dta_chr_sw <- PopGenome::sliding.window.transform(
    dta_chr,
    width = 100,
    jump = 50,
    type = 2
  )

  # Set up sliding windows

  # Chromosome size
  chr_size <- dta_chr@n.sites
  # Set window size and window jump
  window_size <- 100
  window_jump <- 50
  # Use seq to find the start points of each window
  window_start <- seq(from = 1, to = chr_size, by = window_jump)
  # Add the size of the window to each start point
  window_stop <- window_start + window_size
  # Remove windows from the start and stop vectors where the
  # stop position is past the end of the chromosome
  window_start <- window_start[which(window_stop <= chr_size)]
  window_stop <- window_stop[which(window_stop <= chr_size)]
  # Save as a data.frame
  windows_df <- data.frame(
    start = window_start,
    stop = window_stop,
    mid = window_start + (window_stop - window_start) / 2
  )

  # Calculate nucleotide diversity stats (Pi, ...)
  dta_chr_sw <- PopGenome::diversity.stats(dta_chr_sw, pi = TRUE)

  # Calculate FST stats
  dta_chr_sw <- PopGenome::F_ST.stats(dta_chr_sw, mode = "nucleotide")

  # Calculate neutrality stat (Tajima's D, ...)
  dta_chr_sw <- PopGenome::neutrality.stats(dta_chr_sw)

  # Extract statistics for visualization

  # Extract nucleotide diversity, correct for window size
  # and set population names
  nucleotide_diversity <- dta_chr_sw@nuc.diversity.within / window_size
  colnames(nucleotide_diversity) <- paste0(pop_names, "_pi")

  # Extract FST values and set population names
  # Note that here, we need to use t() to transpose the F_ST matrix so that each
  # column is a pairwise comparison and each row is an estimate for a genome
  # window. Since F_ST is pairwise, the column names are also quite different
  # and will also be the same for d_XY_, which is also a pairwise measure.
  fst <- t(dta_chr_sw@nuc.F_ST.pairwise)
  colnames(fst) <- paste0(pop_names[1], "_", pop_names[2], "_fst")

  # Extract dxy (pairwise absolute nucleotide diversity)and set population names
  dxy <- PopGenome::get.diversity(dta_chr_sw, between = TRUE)[[2]] / window_size
  colnames(dxy) <- paste0(pop_names[1], "_", pop_names[2], "_dxy")

  # Combine all statistics into a single data frame
  dta_popgenome_stats <- tibble::as_tibble(data.frame(
    windows_df,
    nucleotide_diversity,
    fst,
    dxy
  ))

  # Visualize the data

  # Select nucleotide diversity data and calculate means
  dta_popgenome_stats |>
    dplyr::select(tidyselect::contains("pi")) |>
    dplyr::summarize_all(c(
      Mean = ~ mean(.),
      Max = ~ max(.),
      Zeros = ~ sum(. == 0)
    ))

  # Gather data for boxplot
  pi_g <- dta_popgenome_stats |>
    dplyr::select(tidyselect::contains("pi")) |>
    tidyr::gather(key = "population", value = "pi")

  # Log transformation
  pi_g$log_pi <- log10(pi_g$pi)

  # Boxplot log-transformed nucleotide diversity (only the values > 0)
  ggplot2::ggplot(
    pi_g,
    ggplot2::aes(population, log_pi, fill = population)
  ) +
    ggplot2::geom_boxplot(color = "black") +
    ggplot2::scale_fill_manual(values = c("red", "blue")) +
    ggplot2::labs(
      x = NULL,
      y = "Log10(pi)",
      title = "Log-transformed nucleotide diversity (with pi > 0 only)",
      subtitle = paste("Chromosome", chrom)
    )

  # Wilcoxon rank-sum test
  wilcox_test_result <- wilcox.test(log_pi ~ population, data = pi_g)
  print(wilcox_test_result)

  # Density plot log-transformed nucleotide diversity
  ggplot2::ggplot(
    pi_g,
    ggplot2::aes(x = log_pi)
  ) +
    ggplot2::geom_density(
      ggplot2::aes(fill = population),
      alpha = 0.5,
      color = "black"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Log10(pi)",
      title = "Log-transformed nucleotide diversity (with pi > 0 only)",
      subtitle = paste("Chromosome", chrom)
    )

  # Visualize FST
  ggplot2::ggplot(
    dta_popgenome_stats,
    ggplot2::aes(mid / 10^6, SWE_ESP_fst)
  ) +
    ggplot2::geom_line(colour = "red") +
    ggplot2::labs(
      x = "Position (Mb)",
      y = expression(italic(F)[ST]),
      title = "FST variation between populations from Spain and Sweden",
      subtitle = paste("Chromosome", chrom)
    )

  # Combine nucleotide diversity, FST and dxy to examine how they
  # co-vary along the genome

  # Select data of interest
  # Fst values smaller than zero are set to zero
  hs_g <- dta_popgenome_stats |>
    dplyr::select(mid, SWE_pi, ESP_pi, SWE_ESP_fst, SWE_ESP_dxy) |>
    dplyr::mutate(dplyr::across(
      c(SWE_ESP_fst, SWE_ESP_dxy),
      ~ ifelse(. < 0, 0, .)
    )) |>
    tidyr::gather(-mid, key = "stat", value = "value")
  # Reorder the levels
  x <- factor(hs_g$stat)
  x <- factor(x, levels(x)[c(3, 1, 4, 2)])
  hs_g$stat <- x

  # Plot with facets
  ggplot2::ggplot(hs_g, ggplot2::aes(mid / 10^6, value, colour = stat)) +
    ggplot2::geom_line() +
    ggplot2::facet_grid(stat ~ ., scales = "free_y") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(
      x = "Position (Mb)",
      title = "FST, Pi + dxy",
      subtitle = paste("Chromosome", chrom)
    )
}
