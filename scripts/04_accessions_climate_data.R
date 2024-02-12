# Load the required packages ----------------------------------------------

# For determining the project root
library(here)
# For transforming, representing, and plotting data
library(tidyverse)
library(ggrepel)
# For printing tables
library(gt)
# For retrieving and plotting the geographical and climatic data
library(geodata)
library(terra)
library(tidyterra)
library(raster)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggspatial)


# Set scientific notation options and default plotting theme --------------

options(scipen = 1)
options(pillar.sigfig = 4)

ggplot2::theme_set(ggplot2::theme_classic())


# Data file paths ---------------------------------------------------------

dta_accessions_file_path <- here::here(
  "data",
  "accessions.tsv")

dta_geo_dir_path <- here::here(
  "data",
  "geodata")

dta_uvb_file_path <- here::here(
  "data",
  "geodata",
  "56459_UVB1_Annual_Mean_UV-B.asc")

dta_pet_06_file_path <- here::here(
  "data",
  "geodata",
  "CHELSA_pet_penman_06_1981-2010_V.2.1.tif")


# Read in the data --------------------------------------------------------

# Read in the accessions data
dta_accessions <- readr::read_tsv(
  file = dta_accessions_file_path,
  col_types = "ffnn")

# Extract a list of the accession coordinates
acc_coords <- dta_accessions[, c("lon", "lat")]

# Read in the bioclimatic data. If the data is not already present at the
# specified path, it will be downloaded.
dta_bio_climate <- geodata::worldclim_global(
  var = "bio",
  res = 0.5,
  path = dta_geo_dir_path)

# Read in the annual mean UVB index data
# Source: https://www.ufz.de/gluv/index.php?en=32367
dta_uvb_idx <- raster::raster(dta_uvb_file_path)

# Read in the potential evapo-transpiration (PET) monthly data (June)
# averaged for 1981-2010
# Source: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2F
dta_pet_06 <- raster::raster(dta_pet_06_file_path)
# Extract the WorldClim bioclimatic data for the accession coordinates
dta_bio_climate_extr <- terra::extract(
  dta_bio_climate,
  acc_coords,
  method = "bilinear",
  df = TRUE
)


# Extract the bioclimatic data --------------------------------------------

# Extract the annual mean UVB index data for the accession coordinates
dta_uvb_idx_extr <- raster::extract(dta_uvb_idx, acc_coords, df = TRUE)

# Extract the potential evapo-transpiration (PET) monthly data (June)
dta_pet_06_extr <- raster::extract(dta_pet_06, acc_coords, df = TRUE)

# Create a tibble of the accessions data and the extracted bioclimatic data
dta_accessions_climate_df <- dplyr::bind_cols(
  dta_accessions |>
    dplyr::rename(
      "Sample" = "sample",
      "Population" = "pop",
      "Latitude" = "lat",
      "Longitude" = "lon"
    ),
  dta_bio_climate_extr |>
    dplyr::select(-1) |>
    dplyr::rename(
      "Annual mean temperature" = "wc2.1_30s_bio_1",
      "Mean diurnal range" = "wc2.1_30s_bio_2",
      "Isothermality" = "wc2.1_30s_bio_3",
      "Temperature seasonality" = "wc2.1_30s_bio_4",
      "Max temperature of warmest month" = "wc2.1_30s_bio_5",
      "Min temperature of coldest month" = "wc2.1_30s_bio_6",
      "Temperature annual range" = "wc2.1_30s_bio_7",
      "Mean temperature of wettest quarter" = "wc2.1_30s_bio_8",
      "Mean temperature of driest quarter" = "wc2.1_30s_bio_9",
      "Mean temperature of warmest quarter" = "wc2.1_30s_bio_10",
      "Mean temperature of coldest quarter" = "wc2.1_30s_bio_11",
      "Annual precipitation" = "wc2.1_30s_bio_12",
      "Precipitation of wettest month" = "wc2.1_30s_bio_13",
      "Precipitation of driest month" = "wc2.1_30s_bio_14",
      "Precipitation seasonality" = "wc2.1_30s_bio_15",
      "Precipitation of wettest quarter" = "wc2.1_30s_bio_16",
      "Precipitation of driest quarter" = "wc2.1_30s_bio_17",
      "Precipitation of warmest quarter" = "wc2.1_30s_bio_18",
      "Precipitation of coldest quarter" = "wc2.1_30s_bio_19"
    ),
  dta_uvb_idx_extr |>
    dplyr::select(-1) |>
    dplyr::rename(
      "UVB index" = "X56459_UVB1_Annual_Mean_UV.B"
    ),
  dta_pet_06_extr |>
    dplyr::select(-1) |>
    dplyr::rename(
      "PET (June)" = "CHELSA_pet_penman_06_1981.2010_V.2.1"
    )
)

# Table output
tbl_accessions_climate <- dta_accessions_climate_df |>
  gt::gt() |>
  gt::opt_stylize() |>
  gt::tab_options(
    table.align = "left",
    table.font.size = "80%"
  ) |>
  gt::fmt_number(
    columns = c(
      `Annual mean temperature`,
      `Mean diurnal range`,
      `Isothermality`,
      `Max temperature of warmest month`,
      `Min temperature of coldest month`,
      `Temperature annual range`,
      `Mean temperature of wettest quarter`,
      `Mean temperature of driest quarter`,
      `Mean temperature of warmest quarter`,
      `Mean temperature of coldest quarter`
    ),
    decimals = 1
  ) |>
  gt::fmt_integer(
    columns = c(
      `Temperature seasonality`,
      `Annual precipitation`,
      `Precipitation of wettest month`,
      `Precipitation of driest month`,
      `Precipitation seasonality`,
      `Precipitation of wettest quarter`,
      `Precipitation of driest quarter`,
      `Precipitation of warmest quarter`,
      `Precipitation of coldest quarter`,
      `UVB index`
    )
  )
tbl_accessions_climate


# Export table tbl-accessions_climate -------------------------------------

# Export the table as a TSV file
readr::write_tsv(
  dta_accessions_climate_df,
  file = here::here("results", "tbl-accessions_climate.tsv"))

# Export the table as a LaTeX file
gt::gtsave(
  tbl_accessions_climate,
  filename = here::here("results", "tbl-accessions_climate.tex")
)


# Create maps of the bioclimatic data at the accession coordinates --------

world_data <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
europe_data <- world_data[world_data$continent == "Europe", ]

# Create a tibble of the points in Europe
europe_points <- dplyr::bind_cols(
  europe_data,
  sf::st_coordinates(sf::st_centroid(sf::st_make_valid(europe_data)$geometry))
)

ggplot2::ggplot(data = europe_data) +
  ggplot2::geom_sf() +
  ggplot2::coord_sf(
    xlim = c(-10, 25),
    ylim = c(35, 70)
  ) +
  ggplot2::geom_point(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2
  ) +
  ggrepel::geom_text_repel(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat, label = sample),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2.5,
    fontface = "bold",
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Spain"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = -2.8,
    nudge_y = -2.2
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Sweden"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = 2.8,
    nudge_y = 4.5
  ) +
  ggspatial::annotation_scale() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      color = gray(0.5),
      linetype = "dashed",
      linewidth = 0.5
    ),
    panel.background = element_rect(fill = "aliceblue")
  ) +
  ggplot2::labs(
    x = "Longitude",
    y = "Latitude",
    title = "Accessions in Spain and Sweden"
  ) +
  ggplot2::theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggplot2::ggplot(data = europe_data) +
  ggplot2::geom_sf() +
  tidyterra::geom_spatraster(
    data = dta_bio_climate,
    mapping = ggplot2::aes(fill = wc2.1_30s_bio_1)
  ) +
  ggplot2::coord_sf(
    xlim = c(-10, 25),
    ylim = c(35, 70)
  ) +
  ggplot2::scale_fill_viridis_c(limits = c(0, 20)) +
  ggplot2::geom_point(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2
  ) +
  ggrepel::geom_text_repel(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat, label = sample),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2.5,
    fontface = "bold",
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Spain"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = -2.8,
    nudge_y = -2.2
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Sweden"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = 2.8,
    nudge_y = 4.5
  ) +
  ggspatial::annotation_scale() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      color = gray(0.5),
      linetype = "dashed",
      linewidth = 0.5
    ),
    panel.background = element_rect(fill = "aliceblue")
  ) +
  ggplot2::labs(
    x = "Longitude",
    y = "Latitude",
    title = "Accessions in Spain and Sweden",
    subtitle = "Annual mean temperature"
  ) +
  ggplot2::theme(
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggplot2::ggplot(data = europe_data) +
  ggplot2::geom_sf() +
  tidyterra::geom_spatraster(
    data = dta_bio_climate,
    mapping = ggplot2::aes(fill = wc2.1_30s_bio_2)
  ) +
  ggplot2::coord_sf(
    xlim = c(-10, 25),
    ylim = c(35, 70)
  ) +
  ggplot2::scale_fill_viridis_c(limits = c(5, 15)) +
  ggplot2::geom_point(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2
  ) +
  ggrepel::geom_text_repel(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat, label = sample),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2.5,
    fontface = "bold",
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Spain"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = -2.8,
    nudge_y = -2.2
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Sweden"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = 2.8,
    nudge_y = 4.5
  ) +
  ggspatial::annotation_scale() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      color = gray(0.5),
      linetype = "dashed",
      linewidth = 0.5
    ),
    panel.background = element_rect(fill = "aliceblue")
  ) +
  ggplot2::labs(
    x = "Longitude",
    y = "Latitude",
    title = "Accessions in Spain and Sweden",
    subtitle = "Mean diurnal range"
  ) +
  ggplot2::theme(
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggplot2::ggplot(data = europe_data) +
  ggplot2::geom_sf() +
  tidyterra::geom_spatraster(
    data = dta_bio_climate,
    mapping = ggplot2::aes(fill = wc2.1_30s_bio_3)
  ) +
  ggplot2::coord_sf(
    xlim = c(-10, 25),
    ylim = c(35, 70)
  ) +
  ggplot2::scale_fill_viridis_c(limits = c(25, 45)) +
  ggplot2::geom_point(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2
  ) +
  ggrepel::geom_text_repel(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat, label = sample),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2.5,
    fontface = "bold",
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Spain"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = -2.8,
    nudge_y = -2.2
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Sweden"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = 2.8,
    nudge_y = 4.5
  ) +
  ggspatial::annotation_scale() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      color = gray(0.5),
      linetype = "dashed",
      linewidth = 0.5
    ),
    panel.background = element_rect(fill = "aliceblue")
  ) +
  ggplot2::labs(
    x = "Longitude",
    y = "Latitude",
    title = "Accessions in Spain and Sweden",
    subtitle = "Isothermality"
  ) +
  ggplot2::theme(
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggplot2::ggplot(data = europe_data) +
  ggplot2::geom_sf() +
  tidyterra::geom_spatraster(
    data = dta_bio_climate,
    mapping = ggplot2::aes(fill = wc2.1_30s_bio_4)
  ) +
  ggplot2::coord_sf(
    xlim = c(-10, 25),
    ylim = c(35, 70)
  ) +
  ggplot2::scale_fill_viridis_c(limits = c(450, 850)) +
  ggplot2::geom_point(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2
  ) +
  ggrepel::geom_text_repel(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat, label = sample),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2.5,
    fontface = "bold",
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Spain"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = -2.8,
    nudge_y = -2.2
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Sweden"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = 2.8,
    nudge_y = 4.5
  ) +
  ggspatial::annotation_scale() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      color = gray(0.5),
      linetype = "dashed",
      linewidth = 0.5
    ),
    panel.background = element_rect(fill = "aliceblue")
  ) +
  ggplot2::labs(
    x = "Longitude",
    y = "Latitude",
    title = "Accessions in Spain and Sweden",
    subtitle = "Temperature seasonality"
  ) +
  ggplot2::theme(
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggplot2::ggplot(data = europe_data) +
  ggplot2::geom_sf() +
  tidyterra::geom_spatraster(
    data = dta_bio_climate,
    mapping = ggplot2::aes(fill = wc2.1_30s_bio_7)
  ) +
  ggplot2::coord_sf(
    xlim = c(-10, 25),
    ylim = c(35, 70)
  ) +
  ggplot2::scale_fill_viridis_c(limits = c(20, 35)) +
  ggplot2::geom_point(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2
  ) +
  ggrepel::geom_text_repel(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat, label = sample),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2.5,
    fontface = "bold",
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Spain"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = -2.8,
    nudge_y = -2.2
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Sweden"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = 2.8,
    nudge_y = 4.5
  ) +
  ggspatial::annotation_scale() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      color = gray(0.5),
      linetype = "dashed",
      linewidth = 0.5
    ),
    panel.background = element_rect(fill = "aliceblue")
  ) +
  ggplot2::labs(
    x = "Longitude",
    y = "Latitude",
    title = "Accessions in Spain and Sweden",
    subtitle = "Temperature annual range"
  ) +
  ggplot2::theme(
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggplot2::ggplot(data = europe_data) +
  ggplot2::geom_sf() +
  tidyterra::geom_spatraster(
    data = dta_bio_climate,
    mapping = ggplot2::aes(fill = wc2.1_30s_bio_12)
  ) +
  ggplot2::coord_sf(
    xlim = c(-10, 25),
    ylim = c(35, 70)
  ) +
  ggplot2::scale_fill_viridis_c(limits = c(300, 1100)) +
  ggplot2::geom_point(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2
  ) +
  ggrepel::geom_text_repel(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat, label = sample),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2.5,
    fontface = "bold",
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Spain"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = -2.8,
    nudge_y = -2.2
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Sweden"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = 2.8,
    nudge_y = 4.5
  ) +
  ggspatial::annotation_scale() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      color = gray(0.5),
      linetype = "dashed",
      linewidth = 0.5
    ),
    panel.background = element_rect(fill = "aliceblue")
  ) +
  ggplot2::labs(
    x = "Longitude",
    y = "Latitude",
    title = "Accessions in Spain and Sweden",
    subtitle = "Annual precipitation"
  ) +
  ggplot2::theme(
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggplot2::ggplot(data = europe_data) +
  ggplot2::geom_sf() +
  tidyterra::geom_spatraster(
    data = dta_bio_climate,
    mapping = ggplot2::aes(fill = wc2.1_30s_bio_15)
  ) +
  ggplot2::coord_sf(
    xlim = c(-10, 25),
    ylim = c(35, 70)
  ) +
  ggplot2::scale_fill_viridis_c(limits = c(0, 100)) +
  ggplot2::geom_point(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2
  ) +
  ggrepel::geom_text_repel(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat, label = sample),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2.5,
    fontface = "bold",
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Spain"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = -2.8,
    nudge_y = -2.2
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Sweden"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = 2.8,
    nudge_y = 4.5
  ) +
  ggspatial::annotation_scale() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      color = gray(0.5),
      linetype = "dashed",
      linewidth = 0.5
    ),
    panel.background = element_rect(fill = "aliceblue")
  ) +
  ggplot2::labs(
    x = "Longitude",
    y = "Latitude",
    title = "Accessions in Spain and Sweden",
    subtitle = "Precipitation seasonality"
  ) +
  ggplot2::theme(
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggplot2::ggplot(data = europe_data) +
  ggplot2::geom_sf() +
  tidyterra::geom_spatraster(
    data = dta_bio_climate,
    mapping = ggplot2::aes(fill = wc2.1_30s_bio_17)
  ) +
  ggplot2::coord_sf(
    xlim = c(-10, 25),
    ylim = c(35, 70)
  ) +
  ggplot2::scale_fill_viridis_c(limits = c(0, 300)) +
  ggplot2::geom_point(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2
  ) +
  ggrepel::geom_text_repel(
    data = dta_accessions,
    mapping = ggplot2::aes(x = lon, y = lat, label = sample),
    color = ifelse(dta_accessions$pop == "ESP", "red", "darkblue"),
    size = 2.5,
    fontface = "bold",
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Spain"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = -2.8,
    nudge_y = -2.2
  ) +
  ggrepel::geom_text_repel(
    data = europe_points |> dplyr::filter(sovereignt == "Sweden"),
    mapping = ggplot2::aes(x = X, y = Y, label = name),
    fontface = "bold",
    min.segment.length = Inf,
    nudge_x = 2.8,
    nudge_y = 4.5
  ) +
  ggspatial::annotation_scale() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      color = gray(0.5),
      linetype = "dashed",
      linewidth = 0.5
    ),
    panel.background = element_rect(fill = "aliceblue")
  ) +
  ggplot2::labs(
    x = "Longitude",
    y = "Latitude",
    title = "Accessions in Spain and Sweden",
    subtitle = "Precipitation of the driest quarter"
  ) +
  ggplot2::theme(
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )
