library(gridExtra)
library(tidyverse)
library(readxl)
library(patchwork)
library(vegan)
library(DSSAT)
library(ggrepel)
library(ggtext)

#' Arachnida traits
traits_arachnida <- read_xlsx('original_data/Arachnida_traits_20240509.xlsx') |>
  mutate(hunters = ifelse(`web` == 1, ifelse(`ambush hunter` == 1, 'web ambush hunter',
                                             ifelse(`active hunter` == 1, 'web active hunter',
                                                    'web')),
                          ifelse(`ambush hunter` == 1, 'ambush hunter',
                                 ifelse(`active hunter` == 1, 'active hunter', 'error')))) |>
  select(species = Species, hunters) |>
  pivot_longer(-1) |>
  unite(form, c('value', 'name'), sep = '_')

#' Arthropoda traits
traits <- read_csv('original_data/Auchenorrhyncha-Heteroptera_traits_20230126.csv') |>
  mutate(form = ifelse(form == 'Calamagrostis', '*Calamagrostis*', form))
#' vegetation
sites <- read_csv('original_data/Experiment_sample-sites_20220208.csv') |>
  left_join(read_csv('enviPCA/evi_pcs.csv')) |>
  mutate(abandoned = 0) |>
  mutate_cond(site %in% c('BAR', 'CER', 'DUB', 'LET', 'POD', 'SVI'), abandoned = 1) #|>

spe <- read_csv('original_data/Experiment_species-data-long_20240514.csv') |>
  mutate(site = ifelse(grepl('_5|6_', id), paste0(site, 2), site)) #|>

spe_wide <- spe |>
  group_by(id, TAXON, site) |>
  summarise(cover = max(cover)) |>
  select(-site) |>
  mutate(cover = sqrt(cover / 100)) |>
  pivot_wider(names_from = TAXON, values_from = cover, values_fill = 0)

first_and_last <- spe |>
  distinct(site, treatment, year) |>
  left_join(sites) |>
  #filter(abandoned == 0) |>
  group_by(site, treatment) |>
  filter(year == min(year) | year == max(year)) |>
  ungroup()

A <- first_and_last |>
  left_join(left_join(spe_wide[1], spe |> distinct(site, treatment, year, id))) |>
  mutate(year = factor(year, levels = 0:4, labels = c(0, 1, 1, 1, 1))) |>
  mutate(treatment = factor(treatment),
         year = as.numeric(year),
         site = gsub('2', '', site))

A_spe <- A |>
  select(id) |>
  left_join(spe_wide) |>
  select(-id)

cap_vse <- capscale(A_spe[-1] ~ treatment:year +
  Condition(site) +
  Condition(year) +
  Condition(treatment),
                    data = A, distance = "bray", sqrt.dist = T)

# -------------------------------------------------------------------------
fit_spe <- envfit(cap_vse,
                  A_spe[colSums(A_spe != 0) > 10],
                  display = 'lc', choices = c(1, 2),
                  correlation = T,
                  perm = 0)

rownames_to_column(data.frame(scores(fit_spe, 'vectors')), 'taxon') |>
  mutate(short_name = paste(str_sub(word(taxon, 1), 1, 3),
                            str_sub(word(taxon, 2), 1, 3))) |>
  mutate(r = fit_spe$vectors$r) |>
  arrange(-r) |>
  slice(1:50) |>
  rename('axis1' = 2, 'axis2' = 3) |>
  mutate_at(c('axis1', 'axis2'), function(x)x * 9) -> fit

guide <- fit |>
  arrange(taxon) |>
  mutate(x = paste0(short_name, ' = ', taxon, ', '))|>
  select(x) |>
  unlist() |>
  paste0(collapse = '')

bind_cols(as_tibble(scores(cap_vse, choices = c(1, 2))$sites), A) |>
  rename('axis1' = 1, 'axis2' = 2) |>
  mutate(sc = paste(axis1, axis2, sep = '_')) |>
  select(sc, site, year, treatment) |>
  arrange(year) |>
  pivot_wider(names_from = year, values_from = sc) |> #print(n =100) |> 
  separate(`1`, into = c('axis1_start', 'axis2_start'), sep = '_') |>
  separate(`2`, into = c('axis1_end', 'axis2_end'), sep = '_') |>
  mutate_at(c('axis1_start', 'axis1_end', 'axis2_start', 'axis2_end'), as.numeric)|>
  mutate(Treatment = factor(treatment, levels = c(0:4),
                            labels = c('Control',
                                       'Mown once (yearly)',
                                       'Mown twice',
                                       '*Rhinanthus* and mown once',
                                       '*Rhinanthus* and mown twice'))) -> meta

expl_by_axes <- ((round(cap_vse$CCA$eig, 2) / cap_vse$tot.chi) * 100)[1:2]

cols1 <- c(`*Calamagrostis*` = 'red',
           Herb = '#CD23F7',
           Graminoid = '#136C2A',
           Predator = 'black')

plants_abbs <- fit |>
  as_tibble() |>
  distinct(taxon, short_name) |>
  mutate(group = 'Vascular plants')

vegetation <- meta |>
  group_by(Treatment) |>
  summarise_if(is.numeric, mean, na.rm = T) |>
  ggplot() +
  geom_segment(data = fit, aes(x = 0, y = 0, xend = axis1, yend = axis2), arrow = arrow(length = unit(.2, 'cm')), col = 'grey86') +
  geom_segment(aes(x = 0, xend = axis1_end, y = 0, yend = axis2_end, col = Treatment),
               size = 1, arrow = arrow()) +
  scale_colour_manual(name = 'Treatment',
                      values = c('black', '#a6bddb', '#0570b0', '#fec44f', '#cc4c02')) +
  ggnewscale::new_scale_color() +
  scale_x_reverse() +
  scale_y_reverse() +
  geom_text_repel(data = fit|>
    mutate(species = taxon) |>
    left_join(traits |>
                mutate(form = factor(form,
                                     levels = names(cols1)))),
                  aes(axis1, axis2, label = short_name,
                      colour = form), max.overlaps = Inf,
                  size = 3) +
  scale_colour_manual(name = 'Life form/Trophic specialization',
                      values = cols1, drop = F) + #, name = 'Growthform/Trophic group') +
  theme_bw() +
  labs(x = 'RDA1',
       y = 'RDA2') +
  theme(legend.position = 'right', #c(1,1),
        panel.grid = element_blank(),
        legend.text = element_markdown(size = 12),
        legend.justification = c(1, 1), legend.background = element_blank(),
        plot.caption = element_textbox_simple(face = 'italic', size = 5))

# -------------------------------------------------------------------------
###########################################################################
#
# Part 2: Arthropoda
#
###########################################################################
# -------------------------------------------------------------------------
sites <- read_csv('original_data\\Experiment_sample-sites_20220208.csv') %>%
  mutate(abandoned = 0) %>%
  mutate_cond(site %in% c('BAR', 'CER', 'DUB', 'LET', 'POD', 'SVI'), abandoned = 1)

f1 <- 'Arthropoda_species-data-combined_2022.xlsx'
f2 <- 'Arthropoda_groups_2022.xlsx'

df_meta <- read_xlsx(paste0('original_data\\', f1), 'spe_Auch_slouceno') %>%
  distinct(plot, site, treatment)
df_step <- df_meta |>
  mutate( #site = ifelse(treatment %in% 5:6, paste0(site, 2), site),
    treatment = ifelse(treatment %in% 5:6, treatment - 2, treatment))


#
# arrange(site, treatment)
# group_by(site) |>
# mutate(kind = factor(n(), levels = 4:5, labels = c('managed', 'abandoned'))) |>
# ungroup() %>%
# select(-source)

dfs <- tibble(group = c('Auchenorchyncha', 'Heteroptera', 'Arachnida'),
              spe = map(list(read_xlsx(paste0('original_data\\', f1), 'spe_Auchenorrhyncha') %>%
                               select(-site, -treatment),
                             read_xlsx(paste0('original_data\\', f1), 'spe_Heteroptera') %>%
                               select(-site, -treatment),
                             read_xlsx(paste0('original_data\\', f1), 'spe_Arachnida') %>%
                               select(-site, -treatment)),
                        function(spe) df_step %>% left_join(spe))) |>
  mutate(spe_no_controls = spe |> map(function(spe) { spe %>% filter(treatment != 0) }))

dfs |>
  pivot_longer(-1) |>
  mutate(model = map(value, ~capscale(sqrt(.x[-c(1:4)]) ~ factor(treatment) + Condition(site),
                                      distance = 'bray',
                                      sqrt.dist = T,
                                      data = .x[1:4]))) |>

  mutate(anova = map2(model, value, ~anova(.x, permutations = how(blocks = .y$site,
                                                                  nperm = 999)))) |>
  filter(name == 'spe') -> dfs2

# -------------------------------------------------------------------------
species_scores <- map2(dfs2$model, dfs2$value,
                       function(.x, .y, rtresh = 0.01, ftresh = 5) {
                         envfit(.x, .y[-c(1:4)],
                                display = 'lc', choices = c(1, 2),
                                correlation = T,
                                perm = 0) -> fit
                         scores(fit, 'vectors') |>
                           as.data.frame() |>
                           rownames_to_column('species') |>
                           as_tibble() |>
                           mutate(r = fit$vectors$r,
                                  freq = colSums(.y[-c(1:4)] != 0)) |>
                           filter(r > rtresh & freq > ftresh) |>
                           rename(axis1 = CAP1, axis2 = CAP2) |>
                           mutate_at(c('axis1', 'axis2'),
                                     function(x) x * 3) |>
                           mutate(short_name = paste(str_sub(word(species, 1), 1, 3),
                                                     str_sub(word(species, 2), 1, 3)))
                       })
tibble(
  abbs = species_scores |>
    append(list(plants_abbs |>
                  rename(species = taxon))) |>
    map_chr(~distinct(.x, species, short_name) |>
      mutate(x = paste0(short_name, ': ', species),
             x = gsub("\\s*\\([^\\)]+\\)", "", as.character(x))) |>
      pull(x) |>
      paste(collapse = ', '))) |>
  write_csv('results/Figure 2 Species Abbreviations.csv')

site_scores <- map2(dfs2$model, dfs2$value,
                    function(.x, .y) {
                      .y <- .y[1:3]
                      bind_cols(as_tibble(scores(.x, choices = c(1, 2))$sites), .y) |>
                        rename('axis1' = 1, 'axis2' = 2) |>
                        mutate(Treatment = factor(treatment, levels = c(0:4),
                                                  labels = c('Control',
                                                             'Mown once (yearly)',
                                                             'Mown twice',
                                                             '*Rhinanthus* and mown once',
                                                             '*Rhinanthus* and mown twice'))) |>
                        group_by(Treatment) |>
                        summarise_at(c('axis1', 'axis2'), mean)
                    })

plots <- map2(site_scores[1:2], species_scores[1:2],
              function(.x, .y) {
                .x |>
                  ggplot(aes(axis1, axis2)) +
                  geom_segment(data = .y, aes(x = 0, y = 0,
                                              xend = axis1,
                                              yend = axis2),
                               arrow = arrow(length = unit(.2, 'cm')), col = 'grey86') +
                  geom_segment(aes(x = 0, xend = axis1,
                                   y = 0, yend = axis2, col = Treatment),
                               size = 1, arrow = arrow()) +
                  scale_colour_manual(values = c('black', '#a6bddb', '#0570b0', '#fec44f', '#cc4c02')) +
                  ggnewscale::new_scale_color() +
                  geom_text_repel(data = .y |>
                    mutate(species) |>
                    left_join(traits |> mutate(form = as.factor(form))),
                                  aes(axis1, axis2, label = short_name,
                                      colour = form), max.overlaps = Inf, size = 3) +
                  scale_colour_manual(values = cols1) +
                  theme_bw() +
                  labs(x = 'RDA1', y = 'RDA2') +
                  theme(legend.position = 'none',
                        legend.text = element_markdown(size = 12),
                        panel.grid = element_blank(),
                        plot.caption = element_textbox_simple(face = 'italic', size = 5))
              })

plots[[3]] <- site_scores[[3]] |>
  ggplot(aes(axis1, axis2)) +
  geom_segment(data = species_scores[[3]], aes(x = 0, y = 0,
                                               xend = axis1,
                                               yend = axis2),
               arrow = arrow(length = unit(.2, 'cm')), col = 'grey86') +
  geom_segment(aes(x = 0, xend = axis1,
                   y = 0, yend = axis2, col = Treatment),
               size = 1, arrow = arrow(),
               show.legend = F) +
  scale_colour_manual(values = c('black', '#a6bddb', '#0570b0', '#fec44f', '#cc4c02')) +
  ggnewscale::new_scale_color() +
  geom_text_repel(data = species_scores[[3]] |>
    mutate(species) |>
    left_join(traits_arachnida |>
                mutate(form =
                         stringr::str_to_title(as.factor(gsub('_hunters', '', form))))),
                  aes(axis1, axis2, label = short_name),
                  max.overlaps = Inf, size = 3) +
  scale_colour_discrete(name = 'Hunting strategies of spiders') +
  #scale_colour_manual(values = cols1) +
  theme_bw() +
  labs(x = 'RDA1', y = 'RDA2') +
  theme( #legend.position = 'none',
    panel.grid = element_blank(),
    plot.caption = element_textbox_simple(face = 'italic', size = 5))

z <- (vegetation + plots[[1]]) / (plots[[2]] + plots[[3]]) +
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = c('a'),
                  tag_suffix = ')',
                  tag_prefix = '(',
                  tag_sep = '') &
  theme(legend.justification = c(1, 1),
        legend.title = element_text(size = 14),
        plot.tag = element_text(size = 18))

z
ggsave('results/Figure 2.svg',
       z,
       width = 13,
       height = 9)
