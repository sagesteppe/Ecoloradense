library(tidyverse)

setwd('~/Documents/Ecoloradense/scripts')

p2tabs <- file.path('..', 'results', 'CountModels', 'tables')
f <- file.path(p2tabs, list.files(p2tabs))

dat <- lapply(f, read.csv)
names(dat) <- gsub('.csv', '', basename(f))
dat <- bind_rows(dat, .id = 'parameters') |>
  mutate(
    Model = str_replace(Model, 'Arithmatic', 'Arithmetic'),
    Resolution = case_when(
      str_detect(parameters, '1-3arc') ~ '1-3arc', 
      str_detect(parameters, '1arc') ~ '1arc', 
      str_detect(parameters, '3arc') ~ '3arc'
    )) |>
  drop_na()

f_lvls <- dat |>
  group_by(Model, Metric) |>
  filter(Metric == 'MAE') |>
  summarize(med = median(Value, na.rm = TRUE)) |>
  arrange(-med) |>
  ungroup() |>
  pull(Model)

dat <- dat |>
  mutate(Model = factor(Model, levels = f_lvls))

rm(p2tabs, f)


# we also will add a 

ggplot(dat, aes(x = Value, y = Model, color = Resolution)) + 
  facet_wrap(~Metric, scales = 'free_x') + 
  geom_jitter(width = 0, height = 0.2) + 
  theme_minimal()

head(dat)

