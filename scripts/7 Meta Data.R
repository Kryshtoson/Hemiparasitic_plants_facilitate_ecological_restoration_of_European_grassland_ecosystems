#' Rhinanthus min mean max cover
read_xlsx('original_data/Calamagrostis-experiment_species-data-long_with_Rhinanthus_20240514.xlsx') |>
  #  filter(grepl('Rhinanthus', TAXON) & grepl('_3_|_4_|_5_|_6_', plot)) |>
  separate(plot, c('site', 'treatment', 'year')) |>
  mutate_at(c('treatment', 'year'), as.numeric) |>
  filter(site != 'ROH') |>
  group_by(site, treatment) |>
  mutate(year = as.numeric(year) - min(as.numeric(year))) |>
  filter(grepl('Rhinanthus', TAXON)) |>
  filter(year > 0) |>
  mutate(treatment = as.numeric(treatment),
         treatment = ifelse(treatment %in% 5:6, treatment - 2, treatment)) |>
  group_by(treatment, year, site) |>
  summarise(cover = sum(cover)) |>
  group_by(treatment, year) |>
  summarise(min = min(cover),
            q05 = quantile(cover, .05),
            q25 = quantile(cover, .25),
            mean = mean(cover),
            q75 = quantile(cover, .75),
            q95 = quantile(cover, .95),
            max = max(cover)) |>
  arrange(year, treatment) |>
  write_xlsx('C:/Users/krystof/OneDrive - MUNI/2022_Calamagrostis/Results/Rhinanthus_cover/Rhinanthus_cover_yearly_summary_stats.xlsx')

read_xlsx('original_data/Calamagrostis-experiment_species-data-long_with_Rhinanthus_20240514.xlsx') |>
  #  filter(grepl('Rhinanthus', TAXON) & grepl('_3_|_4_|_5_|_6_', plot)) |>
  separate(plot, c('site', 'treatment', 'year')) |>
  mutate_at(c('treatment', 'year'), as.numeric) |>
  filter(site != 'ROH') |>
  group_by(site, treatment) |>
  mutate(year = as.numeric(year) - min(as.numeric(year))) |>
  filter(grepl('Rhinanthus', TAXON)) |>
  filter(year > 0) |>
  mutate(treatment = as.numeric(treatment),
         treatment = ifelse(treatment %in% 5:6, treatment - 2, treatment)) |>
  group_by(treatment, year, site) |>
  summarise(cover = sum(cover)) |>
  filter(treatment %in% 3:4) |>
  group_by(site) |>
  summarise(mean = mean(cover)) |>
  write_xlsx('C:/Users/krystof/OneDrive - MUNI/2022_Calamagrostis/Results/Rhinanthus_cover/Rhinanthus_cover_mean_of-3-4.xlsx')