PROJ_ROOT <- '/home/sagesteppe/Documents/Ecoloradense'

library(tidyverse)
library(nlme)
set.seed(27)

oos <- read.csv(file.path(PROJ_ROOT, 'results_summary', 'all_fit_metrics.csv'))
gt_metrics <- read.csv(file.path(PROJ_ROOT, 'results_summary', 'all_fit_2026gt.csv'))

## nAP = (AP - baseline prevalence) / (1 - baseline). OOS (holdout) baseline = the fit's own
## prevalence 1/(1+ratio); GT baseline = the field prevalence of that domain.
nAP <- function(ap, prev) (ap - prev) / (1 - prev)

oos$nAP  <- nAP(oos$pr_auc, 1 / (1 + oos$PAratio))
gt_metrics$nAP <- nAP(gt_metrics$pr_auc, gt_metrics$prevalence)

## One OOS row per model (identity = tree + fname; fname carries res/ratio/seed but NOT the
## tree, so tree is needed to disambiguate designs). Joined onto each of that model's GT rows
## (one per domain x calibration). Keep only nAP and ROC-AUC from each side.
oos_sel <- setNames(oos[, c('tree', 'fname', 'nAP', 'roc_auc')],
                    c('tree', 'fname', 'oos_nAP', 'oos_roc'))

joined <- merge(
  gt_metrics[, c('tree', 'fname', 'resolution', 'PAratio', 'seed', 'domain', 'calibration',
                 'nAP', 'roc_auc')],
  oos_sel, by = c('tree', 'fname'), all.x = TRUE)
names(joined)[names(joined) == 'nAP']     <- 'gt_nAP'
names(joined)[names(joined) == 'roc_auc'] <- 'gt_roc'
## joined: tree, fname, resolution, PAratio, seed, domain, calibration, gt_nAP, gt_roc, oos_nAP, oos_roc


joined = joined |>
 # select(-fname) |>
  filter(domain != 'Extrapolation' & tree != 'Results_spatialknn_iter1a_high_pa' & calibration == 'raw') 


# for spatial mapping grab one seed
test = joined |>
  group_by(tree, resolution, PAratio, domain) |>
  slice_sample(n = 1)

## looks like ecological fallacy. 
nAP = test |>
  select(tree, resolution, PAratio, domain, gt_nAP, oos_nAP) 

model <- lm(gt_nAP ~ oos_nAP + domain + resolution + tree, data = nAP)  # tree = CV strategy
model
summary(model)

## unique contribution of oos_nAP (added last): drop it and compare to the full model
m0 <- update(model, . ~ . - oos_nAP)
partial_r2 <- function(m) { m0 <- update(m, . ~ . - oos_nAP); (deviance(m0) - deviance(m)) / deviance(m0) }
cat(sprintf('oos_nAP: deltaR2 = %.3f, partial R2 = %.3f\n',
            summary(model)$r.squared - summary(m0)$r.squared, partial_r2(model)))

## influence: Cook's distance + leverage. Is the oos_nAP signal driven by a few points?
infl <- transform(nAP, cook = cooks.distance(model), hat = hatvalues(model))
cat(sprintf('Cook cutoff 4/n=%.3f (%d exceed); leverage 2p/n=%.3f (%d exceed)\n',
            4/nrow(nAP), sum(infl$cook > 4/nrow(nAP)),
            2*length(coef(model))/nrow(nAP), sum(infl$hat > 2*length(coef(model))/nrow(nAP))))
print(head(infl[order(-infl$cook), c('tree','resolution','PAratio','gt_nAP','oos_nAP','cook','hat')], 5), row.names = FALSE)
## refit without the top-k most influential -> does oos_nAP slope / partial R2 survive?
for (k in c(2, 5, 10)) {
  mk <- lm(gt_nAP ~ oos_nAP + domain + resolution, data = nAP[-order(-infl$cook)[seq_len(k)], ])
  cat(sprintf('drop top %2d: oos_nAP slope %.3f (full %.3f), partial R2 %.3f (full %.3f)\n',
              k, coef(mk)['oos_nAP'], coef(model)['oos_nAP'], partial_r2(mk), partial_r2(model)))
}


## Results / Discussion sentences, populated live from `model`.
d1   <- drop1(model, test = 'F'); rssf <- deviance(model)
pr2  <- function(t) { ss <- d1[t, 'Sum of Sq']; ss / (rssf + ss) }   # partial R2 of a term
pval <- function(t) d1[t, 'Pr(>F)']
ci   <- confint(model)['oos_nAP', ]
fp   <- function(p) { e <- floor(log10(p)); sprintf('%.1fx10^%d^', p / 10^e, e) }  # a x 10^e
## interaction test: does the holdout->GT slope differ by CV strategy? (n.s. -> report as such)
int_p <- anova(model, update(model, . ~ . + oos_nAP:tree))[2, 'Pr(>F)']

results_txt <- sprintf(
'Using a linear fixed-effects model (n = %d, R^2^ = %.2f, adjusted R^2^ = %.2f), we found strong evidence (partial R^2^ = %.2f, p = %s) of a positive, yet weak, relationship between the nAP of held-out test data and nAP of the independent ground truth data (slope = %.2f, 95%% CI %.2f-%.2f). This model accounted for spatial resolution (partial R^2^ = %.2f, p = %s), prediction domain (partial R^2^ = %.2f, p = %s), and cross-validation strategy (partial R^2^ = %.2f, p = %s). The strength of this relationship did not differ among cross-validation strategies (interaction p = %.2f).',
  nobs(model), summary(model)$r.squared, summary(model)$adj.r.squared, pr2('oos_nAP'), fp(pval('oos_nAP')),
  coef(model)['oos_nAP'], ci[1], ci[2],
  pr2('resolution'), fp(pval('resolution')), pr2('domain'), fp(pval('domain')),
  pr2('tree'), fp(pval('tree')), int_p)

discussion_txt <-
'As expected the relationship between evaluation metrics calculated on held out test data and independent performance are weak, indicating considerable variation in model performance is due to variation in test/train split partitions rather than actual model performance.'

cat('\n## Results\n',     results_txt,    '\n', sep = '')
cat('\n## Discussion\n',  discussion_txt, '\n', sep = '')