library(writexl)
library(nlme)
library(tidyverse)
library(readxl)
library(DSSAT)

#' A comparison of species richness and abundance (number of individuals) per different groups
sites <- read_csv('original_data/Experiment_sample-sites_20220208.csv') %>%
  mutate(abandoned = 0) %>%
  mutate_cond(site %in% c('BAR', 'CER', 'DUB', 'LET', 'POD', 'SVI'), abandoned = 1)
f1 <- 'Arthropoda_species-data-combined_2022.xlsx'
f2 <- 'Arthropoda_groups_2022.xlsx'

bind_rows(
  read_xlsx(paste0('original_data\\', f1), 'spe_Auchenorrhyncha') %>%
    pivot_longer(-(1:3), values_to = 'abundance', names_to = 'species') %>%
    mutate(pa = as.numeric(abundance != 0)) %>%
    mutate(group = 'Auchenorhyncha'),
  read_xlsx(paste0('original_data\\', f1), 'spe_Heteroptera') %>%
    pivot_longer(-(1:3), values_to = 'abundance', names_to = 'species') %>%
    mutate(pa = as.numeric(abundance != 0)) %>%
    mutate(group = 'Heteroptera'),
  read_xlsx(paste0('original_data\\', f1), 'spe_Arachnida') %>%
    pivot_longer(-(1:3), values_to = 'abundance', names_to = 'species') %>%
    mutate(pa = as.numeric(abundance != 0)) %>%
    mutate(group = 'Arachnida')) %>%
  group_by(group, plot, site, treatment) %>%
  summarise(abundance = sum(abundance),
            richness = sum(pa)) %>%
  left_join(sites) %>%
  ungroup() -> input_data

input_data |>
  mutate_at(c('abundance', 'richness'), log1p) %>%
  mutate_cond(treatment == 5, treatment = 3) %>%
  mutate_cond(treatment == 6, treatment = 4) %>%
  mutate(treatment = factor(treatment),
         treatment = relevel(treatment, ref = '1'),
         abandoned = factor(abandoned)) %>%
  group_by(group) %>%
  nest() -> m

m %>% mutate(m_abundance = data %>%
  map(function(df) { lme(abundance ~ treatment,
                         random = ~1 | abandoned / site,
                         data = df) }),
             m_richness = data %>%
               map(function(df) { lme(richness ~ treatment,
                                      random = ~1 | abandoned / site,
                                      data = df) })) -> mods

mods |>
  rename(Abundance = m_abundance, Richness = m_richness) |>
  select(-data) |>
  pivot_longer(-1) |>
  mutate(tidy_df = map(value, broom.mixed::tidy)) |>
  select(group, name, tidy_df)|>
  unnest() |>
  mutate(value = paste0(round(estimate, 4),
                        ', ',
                        ifelse(round(p.value, 3) == 0,
                               '<0.001',
                               round(p.value, 3)))) |>
  select(group, name, term, value) |>
  pivot_wider() |>
  mutate(term = factor(term, levels = c('(Intercept)', paste0('treatment', c(0, 2, 3, 4))),
                       labels = c('Intercept', 'Unmanaged', 'Mown twice', 'Rhinanthus & mown once', 'Rhinanthus & mown twice'))) |>
  write_xlsx('outputs/Coefficients Arthropoda Univariate.xlsx')

mods |>
  rename(Abundance = m_abundance, Richness = m_richness) |>
  select(-data) |>
  pivot_longer(-1) |>
  mutate(anovas = map(value, ~rownames_to_column(as.data.frame(anova(.x)), 'term'))) |>
  select(-value) |>
  unnest() |>
  mutate(value = paste0(
    'F:', '(', numDF, ', ', denDF, ')',
    round(`F-value`, 1),
    ', ',
    ifelse(round(`p-value`, 3) == 0,
           '<0.001',
           round(`p-value`, 3)))) |>
  select(name, term, value)|>
  pivot_wider() |>
  write_xlsx('outputs/Anova Arthropoda Univariate.xlsx')

#' ===========================================================================================
#' plotting Figure 3
#' ===========================================================================================
input_data |>
  mutate_cond(treatment == 5, treatment = 3) %>%
  mutate_cond(treatment == 6, treatment = 4) %>%
  mutate(treatment = factor(treatment),
         abandoned = factor(abandoned)) |>
    pivot_longer(c(abundance, richness)) |>
    mutate(name = factor(name, levels = c('abundance', 'richness'),
                         labels = c('Number of individuals', 'Number of species'))) %>%
  ggplot(aes(treatment, value)) +
  geom_boxplot(
    aes(fill = treatment),
    show.legend = T) +
  scale_y_log10() +
  labs(x = 'Treatment') +
  scale_fill_manual(values = c('grey50', '#a6bddb', '#0570b0', '#fec44f', '#cc4c02'),
                    labels = c('Unmown',
                                       'Mown once (yearly)',
                                       'Mown twice',
                                       '*Rhinanthus* and mown once',
                                       '*Rhinanthus* and mown twice')) +
  facet_grid(name ~ group, scales = 'free_y', switch = 'y') +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0, size = 14),
        strip.text.y = element_text(hjust = .5, size = 14),
        strip.placement = 'outside',
        legend.text = element_markdown(),
        legend.justification = 'left',
        legend.position = 'bottom',
        axis.title.y = element_blank())

ggsave('results/Figure 3.svg', height = 6, width = 8)

