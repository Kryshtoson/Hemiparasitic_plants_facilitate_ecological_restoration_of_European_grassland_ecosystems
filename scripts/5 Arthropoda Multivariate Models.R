library(vegan)
library(tidyverse)
library(readxl)
library(DSSAT)

sites <- read_csv('original_data\\Experiment_sample-sites_20220208.csv') %>%
  mutate(abandoned = 0) %>%
  mutate_cond(site %in% c('BAR', 'CER', 'DUB', 'LET', 'POD', 'SVI'), abandoned = 1)

f1 <- 'Arthropoda_species-data-combined_2022.xlsx'
f2 <- 'Arthropoda_groups_2022.xlsx'

df_meta <- read_xlsx(paste0('original_data\\', f1), 'spe_Auchenorrhyncha') %>%
  distinct(plot, site, treatment)

df_step <- df_meta |>
  mutate(treatment = ifelse(treatment %in% 5:6, treatment - 2, treatment))

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
                                                                  nperm = 999)))) -> z

out <- map(z$anova, ~paste0('F(', paste(.x$Df, collapse = ', '), '): ', round(.x$F[[1]], 2), ', p-value: ',
                            .x$`Pr(>F)`[[1]])) |>
  unlist() |>
  as_tibble() |>
  mutate(name = z$name, group = z$group)

names(z$anova) <- paste0(z$group, '_', z$name)
z