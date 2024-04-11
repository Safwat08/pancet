## Build mega table with all filtered DEGs for Shiny App

# Library
library(Seurat)
library(dplyr)

# Get filtered 
tm <- readRDS('../objects/tabulamuris_ecatlas_pancreas_markers_filtered.rds')
ds <- readRDS('../objects/descartes_ecatlas_pancreas_markers_filtered.rds')
ts <- readRDS('../objects/tabulasapiens_ecatlas_pancreas_markers_filtered.rds')

tm <- tm %>%
  select(-feature_human, -feature) %>%
  mutate(neglogunadjpval = -log(p_val))

max_notinf <- max(tm$neglogunadjpval[!is.infinite(tm$neglogunadjpval)])

tm <- tm %>%
  mutate(neglogunadjpval = ifelse(is.infinite(neglogunadjpval), max_notinf, neglogunadjpval),
         log2FC.neglogunadjpval = neglogunadjpval*avg_log2FC)

tm <- tm %>%
  rename(Features = "feature_updated",
         average.log2fc = "avg_log2FC",
         percent.exp.pancreas.ecs = "pct.1",
         percent.exp.difference = "pct.diff",
         normalized.percent.exp.difference = "norm.pct.diff",
         neg.log.adjusted.p.value = "neglogpval",
         neg.log.unadjustd.p.palue = "neglogunadjpval",
         log2fc.x.neg.log.adjusted.p.value = "log2FC.neglogpval",
         log2fc.x.neg.log.unadjusted.p.value = "log2FC.neglogunadjpval") %>%
  select(Features,
         average.log2fc,
         percent.exp.pancreas.ecs,
         percent.exp.difference,
         normalized.percent.exp.difference,
         neg.log.adjusted.p.value,
         neg.log.unadjustd.p.palue,
         log2fc.x.neg.log.adjusted.p.value,
         log2fc.x.neg.log.unadjusted.p.value)

colnames(tm) <- paste0("TabulaMuris.",colnames(tm))

ds <- ds %>%
  select(-feature) %>%
  mutate(neglogunadjpval = -log(p_val))

max_notinf <- max(tm$neglogunadjpval[!is.infinite(tm$neglogunadjpval)])

ds <- ds %>%
  mutate(neglogunadjpval = ifelse(is.infinite(neglogunadjpval), max_notinf, neglogunadjpval),
         log2FC.neglogunadjpval = neglogunadjpval*avg_log2FC)

ds <- ds %>%
  mutate(neglogunadjpval = ifelse(is.infinite(neglogunadjpval), max_notinf, neglogunadjpval))

ds <- ds %>%
  rename(Features = "feature_updated",
         average.log2fc = "avg_log2FC",
         percent.exp.pancreas.ecs = "pct.1",
         percent.exp.difference = "pct.diff",
         normalized.percent.exp.difference = "norm.pct.diff",
         neg.log.adjusted.p.value = "neglogpval",
         neg.log.unadjustd.p.palue = "neglogunadjpval",
         log2fc.x.neg.log.adjusted.p.value = "log2FC.neglogpval",
         log2fc.x.neg.log.unadjusted.p.value = "log2FC.neglogunadjpval") %>%
  select(Features,
         average.log2fc,
         percent.exp.pancreas.ecs,
         percent.exp.difference,
         normalized.percent.exp.difference,
         neg.log.adjusted.p.value,
         neg.log.unadjustd.p.palue,
         log2fc.x.neg.log.adjusted.p.value,
         log2fc.x.neg.log.unadjusted.p.value)

colnames(ds) <- paste0("Descartes.",colnames(ds))

ts <- ts %>%
  select(-feature) %>%
  mutate(neglogunadjpval = -log(p_val))

max_notinf <- max(ts$neglogunadjpval[!is.infinite(ts$neglogunadjpval)])

ts <- ts %>%
  mutate(neglogunadjpval = ifelse(is.infinite(neglogunadjpval), max_notinf, neglogunadjpval),
         log2FC.neglogunadjpval = neglogunadjpval*avg_log2FC)

ts <- ts %>%
  rename(Features = "feature_updated",
         average.log2fc = "avg_log2FC",
         percent.exp.pancreas.ecs = "pct.1",
         percent.exp.difference = "pct.diff",
         normalized.percent.exp.difference = "norm.pct.diff",
         neg.log.adjusted.p.value = "neglogpval",
         neg.log.unadjustd.p.palue = "neglogunadjpval",
         log2fc.x.neg.log.adjusted.p.value = "log2FC.neglogpval",
         log2fc.x.neg.log.unadjusted.p.value = "log2FC.neglogunadjpval") %>%
  select(Features,
         average.log2fc,
         percent.exp.pancreas.ecs,
         percent.exp.difference,
         normalized.percent.exp.difference,
         neg.log.adjusted.p.value,
         neg.log.unadjustd.p.palue,
         log2fc.x.neg.log.adjusted.p.value,
         log2fc.x.neg.log.unadjusted.p.value)

colnames(ts) <- paste0("TabulaSapiens.",colnames(ts))

data <- full_join(ts,tm, by = join_by(TabulaSapiens.Features == TabulaMuris.Features))
data <- full_join(data,ds, by = join_by(TabulaSapiens.Features == Descartes.Features)) %>%
  mutate(Features = TabulaSapiens.Features) %>% 
  select(-TabulaSapiens.Features)

saveRDS(data, file = '../objects/mega_table.rds')