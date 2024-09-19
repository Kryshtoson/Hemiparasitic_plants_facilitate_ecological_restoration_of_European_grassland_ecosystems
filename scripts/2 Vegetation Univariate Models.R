library(ggtext)
library(sf)
library(tidyverse)
library(broom)
library(nlme)
library(vegan)
library(DSSAT)
library(writexl)
library(broom.mixed)
library(patchwork)
library(readxl)

#' =========================================================
#' Data preparation
#' =========================================================
#' environmental axis
sites <- read_csv('original_data/enviPCA/evi_pcs.csv')

#' species data into wide
spe <- read_csv('original_data\\Experiment_species-data-long_20240514.csv') |>
  mutate(site = ifelse(grepl('_5_|_6_', id), paste0(site, 2), site))

spe_wide <- spe |>
  group_by(id, TAXON, site) |>
  summarise(cover = max(cover)) |>
  select(-site) |>
  mutate(cover = cover / 100) |>
  pivot_wider(names_from = TAXON, values_from = cover, values_fill = 0)

#' biomass, richness, cover of calamagrostis, cover of E1, dissimilarity to target from external source
read_csv('original_data\\Experiment_dry-biomass-weight_20240514.csv') |>
  select(site, block, treatment, year, biomass = DW, community_cover = E1) |>
  full_join(spe |>
              group_by(site, treatment, year) |>
              count(name = 'richness')) |>
  full_join(bind_cols(spe_wide[1], tibble(evenness = diversity(spe_wide[-1]) / log(rowSums(spe_wide[-1] != 0)))) |>
              left_join(spe |> distinct(id, site, treatment, year))) |>
  full_join(spe_wide |>  select(id, 'cala_cover' = "Calamagrostis epigejos"))|>
  left_join(read_csv('meta_data/dissimilarity_to_target.csv')) |>
  select(-id) |>
  mutate(treatment = factor(treatment, levels = c(1, 0, 2, 3, 4))) |>
  left_join(sites) |>
  filter(!(site == 'HIG' & treatment == '0' & year == 3)) |>
  mutate(site = gsub('2', '', site)) -> data


#' =========================================================
#' running all models
#' =========================================================
#' lines left from the manual model selection
anova(lme(log1p(biomass) ~ year * PC1 + treatment * year, random = ~1 | site / treatment, data = data |> filter(!is.na(biomass))))
anova(lme(log1p(cala_cover) ~ year * PC1 + year * PC2 + treatment * year, random = ~1 | site / treatment, data = data))
anova(lme(richness ~ year * PC2 + treatment * year, random = ~1 | site / treatment, data = data))
anova(lme(log(community_cover) ~ year, random = ~1 | site / treatment, data = data))
anova(lme(evenness ~ year * PC1 + year * PC2 + treatment * year, random = ~1 | site / treatment, data = data))
anova(lme(target ~year * PC2 + treatment * year, random = ~1 | site / treatment, data = data |> filter(site != 'KOS')))

term_order <- c('(Intercept)', paste0('treatment', c(0, 2, 3, 4)), 'year', 'PC1', 'PC2',
                paste0('year:', 'treatment', c(0, 2, 3, 4)),
                paste0('year:PC', 1:2))

mods <- tibble(models = list(dw = lme(log1p(biomass) ~ year * PC1 + treatment * year, random = ~1 | site / treatment, data = data |> filter(!is.na(biomass))),
                             cala_cover = lme(log1p(cala_cover) ~ year * PC1 + year * PC2 + treatment * year, random = ~1 | site / treatment, data = data),
                             richness = lme(richness ~  year * PC2 + treatment * year, random = ~1 | site / treatment, data = data),
                             e1 = lme(log(community_cover) ~ year, random = ~1 | site / treatment, data = data),
                             evenness = lme(evenness ~ year * PC1 + year * PC2 + treatment * year, random = ~1 | site / treatment, data = data),
                             target = lme(target ~ year * PC2 + treatment * year, random = ~1 | site / treatment, data = data |> filter(site != 'KOS'))))

mods_full <- tibble(models = list(dw = lme(log1p(biomass) ~ year * PC1 + year * PC2 + treatment * year, random = ~1 | site / treatment, data = data |> filter(!is.na(biomass))),
                             cala_cover = lme(log1p(cala_cover) ~ year * PC1 + year * PC2 + treatment * year, random = ~1 | site / treatment, data = data),
                             richness = lme(richness ~  year * PC1 + year * PC2 + treatment * year, random = ~1 | site / treatment, data = data),
                             e1 = lme(log(community_cover) ~ year * PC1 + year * PC2 + treatment * year, random = ~1 | site / treatment, data = data),
                             evenness = lme(evenness ~ year * PC1 + year * PC2 + treatment * year, random = ~1 | site / treatment, data = data),
                             target = lme(target ~ year * PC1 + year * PC2 + treatment * year, random = ~1 | site / treatment, data = data |> filter(site != 'KOS'))))

mods_full |>
    mutate(
    name = names(models),
    anovas = map(models, ~rownames_to_column(as.data.frame(anova(.x)), 'term'))) |>
  select(-models) |>
  unnest() |>
  mutate(value = paste0(
    'F(', numDF, ', ', denDF, ') = ',
    round(`F-value`, 1),
    '; P ',
    ifelse(round(`p-value`, 3) == 0,
           '< 0.001', paste0('= ', round(`p-value`, 3))))) |>
  select(name, term, value)|>
  pivot_wider() |>
  mutate(term = factor(term, levels = c('(Intercept)', 'year', 'PC1', 'PC2','treatment', 'year:PC1', 'year:PC2', 'year:treatment'))) |>
  arrange(term) |>
  pivot_longer(-1) |>
  arrange(name, term) |>
  relocate(name) |>
  write_xlsx('outputs/Anova Vegetation Univariate Full Models.xlsx')

mods |>
  mutate(
    name = names(models),
    anovas = map(models, ~rownames_to_column(as.data.frame(anova(.x)), 'term'))) |>
  select(-models) |>
  unnest() |>
  mutate(value = paste0(
    'F:', '(', numDF, ', ', denDF, ')',
    round(`F-value`, 1),
    ', ',
    ifelse(round(`p-value`, 3) == 0,
           '<0.001', round(`p-value`, 3)))) |>
  select(name, term, value)|>
  pivot_wider() |>
  mutate(term = factor(term, levels = c('(Intercept)', 'year', 'PC1', 'PC2','treatment', 'year:PC1', 'year:PC2', 'treatment:year'))) |>
  arrange(term) |>
  write_xlsx('outputs/Anova Vegetation Univariate.xlsx')

mods |>
  mutate(tidy_df = map(models, broom::tidy),
         name = names(models)) |>
  select(name, tidy_df)|>
  unnest() |>
  mutate(value = paste0(round(estimate, 4),
                        ', ',
                        ifelse(round(`p.value`, 3) == 0,
           '<0.001', round(`p.value`, 3)))) |>
  select(name, term, value) |>
  pivot_wider() |>
  mutate(term = factor(term, levels = term_order)) |>
  write_xlsx('outputs/Coefficients Vegetation Univariate.xlsx')

#' ===========================================================================================
#' plotting Figure 1
#' ===========================================================================================
#' absolute values
lx <- c("'Dry weight biomass of' ~ italic(Calamagrostis) ~ '(g)'",
        "'Cover of' ~ italic(Calamagrostis) ~ '(%)'",
        "'Number of species'",
        "'Pileou`s evenness'",
        "'Cover of the herb layer (%)'",
        "'Dissimilarity to the target vegetation'")
data |>
  pivot_longer(c(biomass, community_cover, richness, evenness, cala_cover, target)) |>
  select(site, treatment, year, name, value) |>
  group_by(site, name, Treatment = treatment, year) |>
  summarise(mean = mean(value, na.rm = T)) |>
  drop_na() |>
  mutate(Treatment = factor(Treatment, levels = c(0:4),
                            labels = c('Unmown',
                                       'Mown once (yearly)',
                                       'Mown twice',
                                       '*Rhinanthus* and mown once',
                                       '*Rhinanthus* and mown twice'))) |>
  ungroup() |>
  mutate(name = factor(name, levels = c('biomass', 'cala_cover', 'richness', 'evenness', 'community_cover', 'target'),
                       labels = lx)) -> df

df |>
  group_by(name, Treatment, year) |>
  mutate(mean = ifelse(name == 'Cover of Calamagrostis (%)',
                mean*100, mean)) |>
  summarise(
    sem = sd(mean) / sqrt(n()),
    mean = mean(mean)) |>
  ungroup() |>
  ggplot(aes(year, mean)) +
  geom_line(aes(col = Treatment),
            position = position_dodge(width = 0.2), linewidth = 1.2) +
  # scale_y_continuous(labels = scales::percent) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem,
                    col = Treatment), width = 0, linewidth = 1.2,
                position = position_dodge(width = 0.2)) +
  geom_point(aes(col = Treatment),
             size = 3, show.legend = F,
             position = position_dodge(width = 0.2)) +
  scale_colour_manual(values = c('black', '#a6bddb', '#0570b0', '#fec44f', '#cc4c02')) +
  theme_bw() +
  labs(x = 'Experimental year')  +
  facet_wrap(vars(name),
             labeller = label_parsed,
             scales = 'free', strip.position = 'left') +
  theme(legend.justification = c(1, 1),
        strip.placement = 'outside',
        legend.text = element_markdown(),
        strip.background = element_blank(),
        strip.text = element_text(size = 17),
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        text = element_text(size = 15))

ggsave(paste0('results/Figure 1.svg'),
       width = 18,
       height = 10)

#' ==============================================================================
#' relative values
lx <- c("'Dry weight biomass of' ~ italic(Calamagrostis)",
        "'Cover of' ~ italic(Calamagrostis)",
        "'Number of species'",
        "'Pileou`s evenness'",
        "'Cover of the herb layer'",
        "'Dissimilarity to the target vegetation'")
data |>
  pivot_longer(c(biomass, community_cover, richness, evenness, cala_cover, target)) |>
  select(site, treatment, year, name, value) |>
  group_by(site, name, Treatment = treatment, year) |>
  summarise(mean = mean(value, na.rm = T)) |>
  drop_na() |>
  mutate(Treatment = factor(Treatment, levels = c(0:4),
                            labels = c('Unmown',
                                       'Mown once (yearly)',
                                       'Mown twice',
                                       '*Rhinanthus* and mown once',
                                       '*Rhinanthus* and mown twice'))) |>
  ungroup() |>
  mutate(name = factor(name, levels = c('biomass', 'cala_cover', 'richness', 'evenness', 'community_cover', 'target'),
                       labels = lx)) -> df

left_join(df,
          df |>
            filter(year == 0) |>
            select(name, site, Treatment, mean) |>
            rename(mean_y0 = mean)) |>
  mutate(mean_rel = (mean / mean_y0) - 1) |>
  group_by(name, Treatment, year) |>
  summarise(mean_relative = mean(mean_rel),
            sem = sd(mean_rel) / sqrt(n())) |>
  ungroup() |>
  ggplot(aes(year, mean_relative+1)) +
  geom_line(aes(col = Treatment),
            position = position_dodge(width = 0.2), linewidth = 1.2) +
  scale_y_continuous(labels = scales::percent) +
  geom_errorbar(aes(ymin = mean_relative+1 - sem, ymax = mean_relative+1 + sem,
                    col = Treatment), width = 0, linewidth = 1.2,
                position = position_dodge(width = 0.2)) +
  geom_point(aes(col = Treatment),
             size = 3, show.legend = F,
             position = position_dodge(width = 0.2)) +
  scale_colour_manual(values = c('black', '#a6bddb', '#0570b0', '#fec44f', '#cc4c02')) +
  theme_bw() +
  labs(x = 'Experimental year')  +
  facet_wrap(vars(name),
             labeller = label_parsed,
             scales = 'free', strip.position = 'left') +
  theme(legend.justification = c(1, 1),
        legend.text = element_markdown(),
        strip.placement = 'outside',
        strip.background = element_blank(),
        strip.text = element_text(size = 18),
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        text = element_text(size = 15))

ggsave(paste0('results/Figure 1B.svg'),
       width = 18,
       height = 10)