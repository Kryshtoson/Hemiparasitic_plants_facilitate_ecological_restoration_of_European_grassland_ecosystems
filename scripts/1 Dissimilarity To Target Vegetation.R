library(sf)
library(tidyverse)
library(broom)
library(nlme)
library(vegan)
library(DSSAT)
library(writexl)
library(broom.mixed)
library(readxl)
library(patchwork)

#' harmonized species list with red-list status for Czech Republic and Slovakia separately and header data
species_list <- read_csv('original_data\\Target-vegetation_species-list_20220223.csv') %>%
  select(species, new_species, keep, Cesko, Slovensko)
head <- read_xlsx('original_data\\Target-vegetation_20240517.xlsx')
country_id <- head |>
  select(reference = id, name = Source) |>
  filter(name %in% c('CZ', 'SK'))

#read_csv('original_data\\Target-vegetation_species-data-long_20220223.csv')

#' species data from Czech and Slovak national dabases (background data)
spe_long <- read_csv('original_data\\Target-vegetation_species-data-long_20220223.csv') %>%
  left_join(species_list) %>%
  filter(keep) %>%
  mutate(species = new_species) %>%
  group_by(id = `Releve number`, species) %>%
  summarise(cover = sum(cover))

spe_wide <- spe_long |>
  pivot_wider(names_from = species, values_from = cover, values_fill = 0) |>
  column_to_rownames('id') |>
  mutate_all(sqrt)

# Dissimilarity Calculation:
#
# The species data is first transformed from long to wide format,
# and then a dissimilarity matrix is calculated using the Bray-Curtis method.
# A distance matrix based on geographic coordinates is also calculated after
# converting the coordinates to a metric system using the sf package.
# This involves spatial transformations and the use of st_distance() for calculating distances.

# Only plots within a 5000 m radius are selected for comparison, except for Slovakian plots,
# where insufficient number of target plots made us combine all target plots for all experimental sites.
# This is done by mannually adjusting the spatial distances between plots.

species_list |>
  filter(keep) |>
  select(species = new_species, CZ = Cesko, SVK = Slovensko) |>
  pivot_longer(c(CZ, SVK)) |>
  count(value)

no_redlisted <- spe_long |>
  rename(reference = id) |>
  left_join(country_id) |>
  mutate(reference = as.character(reference)) |>
  left_join(
    species_list |>
      filter(keep) |>
      select(species = new_species, CZ = Cesko, SVK = Slovensko) |>
      pivot_longer(c(CZ, SVK)) |>
      filter(value == c('CR', 'EN', 'NT', 'VU'))) |>
  summarise(redlisted = sum(!is.na(value)))
distm_spe <- vegdist(spe_wide, 'bray')

spatial_dist <- head |>
  st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326) |>
  st_transform(crs = 25833) |>
  mutate(X = st_coordinates(geometry)[, 'X'],
         Y = st_coordinates(geometry)[, 'Y']) |>
  as_tibble() |>
  column_to_rownames('id') |>
  select(X, Y) |>
  dist() |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column('id') |>
  as_tibble() |>
  pivot_longer(-id, names_to = 'reference', values_to = 'distance') |>
  filter(id != reference) |>
  mutate(id = as.numeric(id),
         distance = ifelse(reference %in% as.character(country_id$reference[country_id$name == 'SK']) &
                             id %in% head$id[head$Site %in% c('ILI', 'KYS', 'DUC')], 99, distance)) |>
  filter(distance < 10000)

compositional_dist <- distm_spe |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column('id') |>
  as_tibble() |>
  pivot_longer(-id, names_to = 'reference', values_to = 'dissimilarity') |>
  filter(id != reference) |>
  mutate(id = as.numeric(id)) |>
  semi_join(head |> filter(Source == 'Orig')) |>
  filter(!reference %in% id) |>
  semi_join(spatial_dist) |>
  left_join(no_redlisted) |>
  group_by(id) |>
  arrange(-redlisted) |>
  slice(1:(n() * .6))

# Number of reference plots per site
compositional_dist |>
  left_join(head) |>
  group_by(Site, id) |>
  count() |>
  group_by(Site) |>
  summarise(n = mean(n)) |>
  print(n = 100)

compositional_dist |>
  group_by(id) |>
  summarise(target = min(dissimilarity)) |>
  left_join(head[1:2]) |>
  select(id = Reference, target) |>
  write_csv('meta_data/dissimilarity_to_target.csv')

# number of plots in the dataset from Czech Republic and Slovakia
compositional_dist |>
  mutate(reference = as.numeric(reference)) |>
  left_join(country_id) |>
  ungroup() |>
  distinct(name, reference) |>
  count(name)
