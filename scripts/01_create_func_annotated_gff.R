# Load libraries ----------------------------------------------------------

# For determining the project root
library(here)
# For transforming and representing data
library(tidyverse)

# Determine the file paths of the data relative to the project's base dir -

dta_gff_file_path <- here::here(
  "data",
  "TAIR10_GFF3_genes.gff"
)

dta_gff_annot_file_path <- here::here(
  "data",
  "TAIR10_GFF3_genes_annot.gff"
)

dta_func_descr_file_path <- here::here(
  "data",
  "TAIR10_functional_descriptions.txt"
)


# Read in the data --------------------------------------------------------

dta_gff <- readr::read_tsv(
  file = dta_gff_file_path,
  col_names = FALSE
)

dta_func_descr <- readr::read_tsv(
  file = dta_func_descr_file_path
)


# Annotate the GFF with the functional annotation--------------------------

# Filter the GFF for mRNA entries and the chromosomes 1 to 5 +
# Extract the chromosome number1 1 to 5 (from Chr1, Chr2, ...) +
# Add Model_name column +
# Join with the functional description file via column Model_name +
# Export the annotated GFF
dta_gff |>
  dplyr::filter(
    X1 %in% c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5") &
      X3 == "mRNA"
  ) |>
  dplyr::mutate(X1 = stringr::str_extract(
    string = X1,
    pattern = "\\d+"
  )) |>
  dplyr::mutate(Model_name = stringr::str_extract(
    string = X9,
    pattern = ".*Name=([\\w\\.]+).*",
    group = 1
  )) |>
  dplyr::left_join(dta_func_descr, by = "Model_name") |>
  readr::write_tsv(file = dta_gff_annot_file_path, col_names = FALSE)
