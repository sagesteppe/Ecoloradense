library(tidyverse)

setwd('~/Documents/Ecoloradense/scripts')

p2tabs <- file.path('..', 'results', 'CountModels', 'tables')
f <- file.path(p2tabs, list.files(p2tabs))
cols <- c("3arc" = "#0FA3B1", "1arc" = "#b07e62", "1-3arc" = "#40501e")
labs <- c("3arc" = "3 arc", "1arc" = "1 arc", "1-3arc" = "1/3 arc")


suit <- lapply(f, read.csv)
names(suit) <- gsub('.csv', '', basename(f))
suit <- bind_rows(suit, .id = 'parameters') |>
  mutate(
    Model = str_replace(Model, 'Arithmatic', 'Arithmetic'),
    Resolution = case_when(
      str_detect(parameters, '1-3arc') ~ '1-3arc', 
      str_detect(parameters, '1arc') ~ '1arc', 
      str_detect(parameters, '3arc') ~ '3arc'
    )) |>
  drop_na()

f_lvls <- suit |>
  group_by(Model, Metric) |>
  filter(Metric == 'MAE') |>
  summarize(med = median(Value, na.rm = TRUE)) |>
  arrange(-med) |>
  ungroup() |>
  pull(Model)

suit <- suit |>
  mutate(Model = factor(Model, levels = f_lvls))

rm(p2tabs, f)

# we also will add a 

ggplot(suit, aes(x = Value, y = Model, color = Resolution)) + 
  facet_wrap(~Metric, scales = 'free_x') + 
  scale_color_manual(values = cols, labels = labs) + 
  geom_jitter(width = 0, height = 0.2) + 
  theme_minimal() +  
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5)) + 
  theme(
    legend.position = 'bottom',
    legend.background = element_rect(color = "grey75"),
    plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'Evaluation metrics for count prediction models')

ggsave( file.path('..', 'figures', 'CountPredictions.png'), height = 6, width = 6)

#####################################################

# take a look at the prediction data. 

p2tabs <- file.path('..', 'results', 'tables')
f <- file.path(p2tabs, list.files(p2tabs, pattern = 'csv'))

auc <- lapply(f, read.csv)
names(auc) <- gsub('.csv', '', basename(f))

auc <- bind_rows(auc, .id = 'run') 

f1 <- file.path(p2tabs, list.files(p2tabs, pattern = 'rds'))
ksink <- lapply(f1, readRDS)
names(ksink) <- gsub('.csv', '', basename(f1))

ksink <- left_join(
  lapply(ksink, '[[', 'overall') |>
    bind_rows(.id = 'run'), 
  lapply(ksink, '[[', 'byClass') |>
    bind_rows(.id = 'run'), 
  by = 'run'
) |>
  pivot_longer(!run, values_to = 'estimate', names_to = 'metric') |>
  bind_rows(select(auc, run, metric, estimate)) 

rm(f, f1, auc, p2tabs)


# clean up data for plotting

ksink <- ksink |>
  mutate(
    resolution = str_extract(run, '.*arc'), 
    PA_ratio = as.numeric(str_remove_all(str_extract(run, ':.*D'), ':|D')),
    metric = str_replace(metric, 'pr_auc', 'PR AUC'), 
    metric = str_replace(metric, 'roc_auc', 'ROC AUC'), 
    metric = str_replace(metric, 'Neg Pred Value', 'Negative Pred. Value'),
    metric = str_replace(metric, 'Neg Pred Value', 'Positive Pred. Value'),
  ) |>
  select(-run) |>
  filter(! metric %in% c(
    'AccuracyLower', 'AccuracyUpper',
    'AccuracyNull', 'AccuracyPValue', 'McnemarPValue')
    ) 

table(ksink$metric)

sink_lvls <- c(
  'Balanced Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'Recall', 
  'Kappa', 'PR AUC', 'ROC AUC', 'F1')

ksink |>
  mutate(metric = forcats::fct_relevel(metric, sink_lvls)) |>
  ggplot( aes(x = estimate, y = resolution, color = resolution)) + 
  facet_wrap(~metric, scales = 'free_x', nrow = 5) + 
  scale_color_manual(values = cols, labels = labs, name = 'Resolution') + 
  geom_point() + 
  theme_bw() + 
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5)) + 
  theme(
    legend.position = 'bottom',
    legend.background = element_rect(color = "grey75"),
    plot.title = element_text(hjust = 0.5)) + 
  labs(
    title = 'Prediction of the probability of suitable habitat', 
    x = 'Estimate', y = 'Resolution')


ggsave( file.path('..', 'figures', 'PrSuitPredictions.png'), height = 9, width = 6)

# paradoxical results - may be due to overall decrease in independent variable
# accuracy from native grain (1arc decrease), while 1-3 arc overcomes the limitation by 
# sample size? ?? 

# may need to utilize the same indx data splits.... ?? 


# in terms of PA ratios only trends seem to be that higher ratio is bad for 1-3, 
# and this might be flipped for 3 arc. 

# could be worth running a 2.6 and a 3.6 to increase this range and see if a trend exists. 