# Load required packages
library(patchwork)
library(ggnewscale)
library(tidyverse)
library(vegan)
library(knitr)
library(DSSAT)
library(ggrepel)

#' =====================================================================================
#' Data input
#' Assignment of abandoned plots
#' Conversion of plots to wide form
#' ballancing Rhinanthus minor treatment (coded as 5 and 6) by replicating treatments 1:3 on ILI, KYS and DUC
sites <- read_csv('original_data/Experiment_sample-sites_20220208.csv') |>
  left_join(read_csv('original_data/enviPCA/evi_pcs.csv')) |>
  mutate(abandoned = 0) |>
  mutate_cond(site %in% c('BAR', 'CER', 'DUB', 'LET', 'POD', 'SVI'), abandoned = 1)

spe <- read_csv('original_data/Experiment_species-data-long_20240514.csv') |>
  mutate(site = ifelse(grepl('_5|6_', id), paste0(site, 2), site))

spe_wide <- spe |>
  group_by(id, TAXON, site) |>
  summarise(cover = max(cover)) |>
  select(-site) |>
  mutate(cover = sqrt(cover / 100)) |>
  pivot_wider(names_from = TAXON, values_from = cover, values_fill = 0)

first_and_last <- spe |>
  distinct(site, treatment, year) |>
  left_join(sites) |>
  group_by(site, treatment) |>
  filter(year == min(year) | year == max(year)) |>
  ungroup()

#' =====================================================================================
#' Data subsets
#' Analyzed plots are only of the first and the last year of the experiment
#' Defining A, B and C groups
#' A includes all plots just as we got them
#' B does not include control plots (treatment == 0)
#' C includes only those sites where all treatments are available
A_meta <- first_and_last |>
  left_join(left_join(spe_wide[1], spe |> distinct(site, treatment, year, id))) |>
  mutate(year = factor(year, levels = 0:4, labels = c(0, 1, 1, 1, 1))) |>
  mutate(treatment = factor(treatment),
         year = as.numeric(year))
         #site = gsub('2', '', site))
A <- A_meta |>
  bind_rows(A_meta[A_meta$treatment %in% 0:2 & A_meta$site %in% c('DUC', 'ILI', 'KYS'),] |>
  mutate(site = paste0(site, 2)))

B <- A |>
  filter(treatment != 0)

C <- A |>
  filter(!abandoned)

#' selection of species data to be also done for every group separately
A_spe <- A |>
  select(id) |>
  left_join(spe_wide) |>
  select(-id)
B_spe <- B |>
  select(id) |>
  left_join(spe_wide) |>
  select(-id)
C_spe <- C |>
  select(id) |>
  left_join(spe_wide) |>
  select(-id)

#' ===========================================================================
#' modelling
#' ===========================================================================
#' Model A with all original data in
A |>
  count(treatment, site) |>
  ggplot(aes(site, treatment)) +
  geom_label(aes(label = paste('n =', n/2), fill = factor(n)),
             show.legend = F) +
  coord_flip()

A_cap <- capscale(A_spe ~ treatment:year +
  Condition(site + treatment + year * PC1 + year * PC2),
                  data = A,
                  distance = "bray",
                  sqrt.dist = T)

anova(A_cap, permutations = how(blocks = A$site,
                                plots = Plots(paste(A$site, A$treatment),
                                              type = 'free'),
                                within = Within(type = "none"),
                                nperm = 999))

#' Model B with balanced design. Controls removed.
B |>
  count(treatment, site) |>
  ggplot(aes(site, treatment)) +
  geom_label(aes(label = paste('n =', n/2), fill = factor(n)),
             show.legend = F) +
  coord_flip()

B_cap <- capscale(B_spe ~ treatment:year + Condition(site + treatment + year * PC1 + year * PC2),
                  data = B, distance = "bray", sqrt.dist = T)

anova(B_cap, permutations = how(blocks = B$site,
                                plots = Plots(paste(B$site, B$treatment),
                                              type = 'free'),
                                within = Within(type = "none"),
                                nperm = 999))

#' Model C with balanced design. Only sites with all treatments.
C |>
  count(treatment, site) |>
  ggplot(aes(site, treatment)) +
  geom_label(aes(label = paste('n =', n/2), fill = factor(n)),
             show.legend = F) +
  coord_flip()

C_cap <- capscale(C_spe ~ treatment:year + Condition(site + treatment + year * PC1 + year * PC2),
                  data = C, distance = "bray", sqrt.dist = T)

anova(C_cap, permutations = how(blocks = C$site,
                                plots = Plots(paste(C$site, C$treatment),
                                              type = 'free'),
                                within = Within(type = "none"),
                                nperm = 999))