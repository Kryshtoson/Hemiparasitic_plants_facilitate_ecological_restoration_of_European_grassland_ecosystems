library(ggtext)
library(writexl)
library(nlme)
library(readxl)
library(tidyverse)
library(sf)

read_xlsx('original_data/Arachnida_traits_20240509.xlsx') |>
  mutate(structure = ifelse(ground == 1, ifelse(vegetation == 1, 'both', 'ground'),
                            'vegetation')) |>
  mutate(hunters = ifelse(`web` == 1, ifelse(`ambush hunter` == 1, 'web ambush hunter',
                                             ifelse(`active hunter` == 1, 'web active hunter',
                                                    'web')),
                          ifelse(`ambush hunter` == 1, 'ambush hunter',
                                 ifelse(`active hunter` == 1, 'active hunter', 'error')))) |>
  select(species = Species, structure, hunters) |>
  pivot_longer(-1) |>
  unite(form, c('value', 'name'), sep = '_') |>
  bind_rows(read_csv('original_data/Auchenorrhyncha-Heteroptera_traits_20230126.csv')) |>
  group_by(species) |>
  mutate(n = n()) |>
  filter(!(n == 3 & form == 'Predator')) |>
  select(-n) -> arachnida_traits

read_csv('original_data/experiment_species-data-long_20240514.csv') |>
  left_join(read_csv('original_data/Auchenorrhyncha-Heteroptera_traits_20230126.csv'),
            by = c(TAXON = 'species')) |>
  mutate(
    site = factor(ifelse(grepl('_5_|_6_', id), paste(site, 2), site)),
    #treatment = factor(ifelse(grepl('_5_|_6_', id), treatment - 2, treatment)),
    year = factor(year, levels = 0:3),
    treatment = factor(treatment, levels = c(1, 0, 2, 3, 4)),
    form = fct_relevel(factor(form), 'Graminoid')) |>
  group_by(year = as.numeric(year), site, treatment, form) |>
  summarise(cover_sum = sum(cover), drop = F) |>
  group_by(year, site, treatment) |>
  mutate(proportion = cover_sum / sum(cover_sum)) |>
  filter(form != 'Calamagrostis') |>
  group_by(site) |>
  mutate(plot = paste0(site, year, treatment))  |>
ungroup() -> guildy_kytky

m_kytky <- lme(log(cover_sum) ~ form * treatment * year, random = ~1 | site / treatment / plot,
               data = guildy_kytky |> mutate(site = str_sub(site, 1, 3)))
anova(m_kytky)

broom.mixed::tidy(m_kytky) |>
  write_xlsx('outputs/Coefficients Vegetation Guilds.xlsx')

f1 <- 'Arthropoda_species-data-combined_2022.xlsx'
tibble(group = c('Auchenorrhyncha', 'Heteroptera', 'Arachnida'),
       spe = list(read_xlsx(paste0('original_data\\', f1), 'spe_Auchenorrhyncha'),
                  read_xlsx(paste0('original_data\\', f1), 'spe_Heteroptera'),
                  read_xlsx(paste0('original_data\\', f1), 'spe_Arachnida'))) |>
  unnest(spe) |>
  pivot_longer(-c(group, plot, site, treatment)) |>
  mutate(
    site = ifelse(treatment %in% 5:6, paste(site, '2'), site),
    treatment = ifelse(treatment %in% 5:6, treatment - 2, treatment)) |>
  drop_na() |>
  filter(value != 0) |>
  select(-plot) |>
  left_join(arachnida_traits, by = c(name = 'species')) |>
  filter(!is.na(form)) |>
  separate(form, c('form', 'Ara_type'), sep = '_') |>
  mutate(group = ifelse(!is.na(Ara_type), paste(group, Ara_type), group)) |>
  mutate(treatment = factor(treatment, levels = c(1, 0, 2, 3, 4)),
       form = fct_relevel(factor(form), 'Graminoid'),
       group = factor(group),
       site = factor(site)) |>
  group_by(group, treatment, form, site) |>
  summarise(cover_sum = sum(value)) |>
  arrange(group, form) |>
  group_by(group, site, treatment) |>
  mutate(proportion = cover_sum / sum(cover_sum)) -> guildy_hmyz

auch <- guildy_hmyz |> filter(group == 'Auchenorrhyncha')
hete <- guildy_hmyz |>
  filter(group == 'Heteroptera') |>
  filter(!(form == 'Conifer' | form == 'Moss'))
ara_hun <- guildy_hmyz |>
  filter(group == 'Arachnida hunters')
ara_stru <- guildy_hmyz |>
  filter(group == 'Arachnida structure')

table(guildy_hmyz$group)
m_auch <- lme(log1p(cover_sum) ~ form * treatment, random = ~1 | site / treatment, data = auch |> mutate(site = str_sub(site, 1, 3)))
m_hete <- lme(log1p(cover_sum) ~ form * treatment, random = ~1 | site / treatment, data = hete |> mutate(site = str_sub(site, 1, 3)))
m_arahu <- lme(log1p(cover_sum) ~ form * treatment, random = ~1 | site / treatment, data = ara_hun |> mutate(site = str_sub(site, 1, 3)))
m_arastr <- lme(log1p(cover_sum) ~ form * treatment, random = ~1 | site / treatment, data = ara_stru |> mutate(site = str_sub(site, 1, 3)))

summary(m_auch)
summary(m_hete)

anova(m_auch)
anova(m_hete)
anova(m_arahu)
anova(m_arastr)

tibble(
  group = c('Auchenorhyncha', 'Heteroptera', 'Arachnida'),
  models = list(m_auch, m_hete, m_arahu)) |>
  mutate(
    name = names(models),
    anovas = map(models, ~rownames_to_column(as.data.frame(anova(.x)), 'term'))) |>
  select(-models) |>
  unnest() |>
  mutate(value = paste0(
    'F:', '(', numDF, ', ', denDF, ')',
    round(`F-value`, 1),
    ', ',
    round(`p-value`, 4))) |>
  select(name = group, term, value) |>
  pivot_wider() |>
  write_xlsx('outputs/Anova Arthropoda Guilds.xlsx')


broom.mixed::tidy(m_hete) |>
  mutate(p.value < 0.05) |>
  select(-effect, -std.error, -statistic)

broom.mixed::tidy(m_arahu) |>
  mutate(p.value = p.value < 0.05)

list(broom.mixed::tidy(m_auch),
     broom.mixed::tidy(m_hete),
     broom.mixed::tidy(m_arahu),
     broom.mixed::tidy(m_arastr)) |>
  setNames(c('Auchenorhyncha', 'Heteroptera', 'Arachnida_hunters', 'Arachnida_structure')) |>
  write_xlsx('outputs/Coefficients Arthropoda Guilds.xlsx')


cols1 <- c(`*Calamagrostis*` = '#f97971',
           Herb = '#9091c8',
           Graminoid = '#A1B65E',
           Conifer = '#9CCB69',
           Moss = '#136C2A',
           Predator = 'grey40')

guildy_hmyz |>
  filter(group != 'Arachnida structure' & group != 'Arachnida hunters') |>
  filter(!(form == 'Conifer' | form == 'Moss')) |>
  mutate(Treatment = factor(treatment, levels = c(0:4),
                            labels = c('Unmown',
                                       'Mown once (yearly)',
                                       'Mown twice',
                                       '*Rhinanthus* and mown once',
                                       '*Rhinanthus* and mown twice')),
         form = factor(form, levels = c('Calamagrostis', 'Graminoid', 'Herb', 'Predator'),
                       labels = c('*Calamagrostis*', 'Graminoid', 'Herb', 'Predator'))) |>
  group_by(group, form, Treatment) |>
  ggplot(aes(Treatment, cover_sum)) +
  geom_boxplot(aes(fill = form), alpha = 1) +
  coord_flip() +
  facet_wrap(~group) +
  scale_fill_manual(name = 'Trophic specialization',
                    values = cols1, drop = F) +
  scale_colour_manual(name = 'Trophic specialization',
                      values = cols1, drop = F) +
  scale_y_log10(expand = c(0, 0, .02, .02)) +
  scale_x_discrete(limits = rev) +
  labs(y = 'Number of individuals') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.justification = c(1, 1),
        strip.placement = 'outside',
        legend.text = element_markdown(),
        axis.text.y = element_markdown(),
        strip.background = element_blank(),
        strip.text = element_text(size = 18, hjust = 0),
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        text = element_text(size = 15))

ggsave('results/Figure 4.svg',
       width = 12, height = 7)

