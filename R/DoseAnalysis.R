# Dose Response Data Analysis
# All wells contain the total binding (TOTB) control which generates the maximum raw data signal, but represents vehicle (0% biological activity).

#Set up -------------------------------------------------
library(tidyverse)
library(lubridate)
library(viridis)
library(ggpubr)
library(rstatix)
library(ROCR)
library(ggplot2)
library(patchwork)

# Prepare data -------------------
DRPlateMap <- read_csv('Data/Platemap.csv')

PlateId2  <- data.frame(PlateId2 = 1:28)

DoseData <- read_csv('Data/AllData.csv')%>%
  filter(PlateMap == 'Dose') %>%
  select(-PlateMap) %>%
  mutate(Assay = factor(Assay, levels = c("Tgt1", "Tgt2")),
         PlateType = as.character(PlateType),
         Well = as_factor(Well),
         PlateId = as.character(PlateId)) %>%
  left_join(DRPlateMap) %>%
  mutate(Sample = fct_cross(Cmpd, as_factor(Conc))) %>%
  group_by(PlateId,ExpTime, Assay, PlateType) %>%
  nest() %>%
  mutate(ExpTime = mdy_hms(ExpTime),
         ExpDate = as_date(ExpTime),
         Run = as_factor(if_else(ExpDate < mdy('12/1/2019'), 'Run1', 'Run2'))) %>%
  arrange(Assay, Run) %>%
  bind_cols(PlateId2)%>%
  mutate(AssayPlate = if_else(Assay == 'Tgt1', PlateId2 - 14, PlateId2 - 0),
         AssayPlate = as_factor(AssayPlate))%>%
  select(-PlateId2) %>%
  unnest(cols= c(data)) %>%
  ungroup()

#  Normalize to Plate controls -------------------
CtrlData <- DoseData %>%
  filter(Cmpd == 'TOTB' | Cmpd == 'NSB') %>%
  group_by(Assay, PlateId, Cmpd) %>%
  summarise(Mean = mean(Data),
            Median = median(Data)) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with('M'), names_to = 'NormStatType', values_to = 'StatVal') %>%
  mutate(NormStat = paste(Cmpd, NormStatType, sep = '_')) %>%
  select(-Cmpd, -NormStatType) %>%
  pivot_wider(names_from = starts_with('Norm'), values_from = StatVal)

DoseData <- DoseData %>%
  left_join(CtrlData) %>%
  mutate(PctActMean = 100 * (1 - ((NSB_Mean - Data) / (NSB_Mean - TOTB_Mean))),
         PctActMedian = 100 * (1 - ((NSB_Median - Data) / (NSB_Median - TOTB_Median)))
  ) %>%
  select(-starts_with('NSB'), -starts_with('TOTB'))

# Plate Controls and QC ---------------------------
PlateQCData <- DoseData %>%
  filter(Cmpd == 'TOTB' | Cmpd == 'NSB') %>%
  select(-Conc, -Sample) %>%
  pivot_longer(cols = starts_with('Pct'), names_to = 'Scale', values_to = 'Activity') %>%
  mutate(Scale = if_else(str_detect(Scale, 'Mean'), 'Mean', 'Median')
         )

# Create Control Well Plots Fig 3
T1Mean <- ggplot(filter(PlateQCData, Assay == 'Tgt1', Scale == 'Mean'), aes(x = PlateId, y = Activity, fill = Cmpd)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, alpha = 0.5) +
  labs(x = 'Plate',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = 'none') +
  labs(title = 'Tgt1')

T2Mean <- ggplot(filter(PlateQCData, Assay == 'Tgt2', Scale == 'Mean'), aes(x = PlateId, y = Activity, fill = Cmpd)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, alpha = 0.5) +
  labs(x = 'Plate',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = 'Tgt2')

Fig3 <- (T1Mean / T2Mean) +
  plot_annotation(title = 'Figure 3 Plate Controls.', tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave('Figures/Weidner Fig 3.jpg', plot = Fig3, height = 6, width = 6, units = 'in', dpi = 300)

# Comparison of Standard and Robust Z-factors Fig 4
PlateQCData <- PlateQCData %>%
  group_by(Assay, Run, AssayPlate, Cmpd, Scale) %>%
  summarise(Avg = mean(Activity),
            Med = median(Activity),
            StDev = sd(Activity),
            MAD = mad(Activity)) %>%
  ungroup() %>%
  mutate(ZLim = if_else(Cmpd == 'NSB', Avg - 3 * StDev, Avg + 3 * StDev),
         ZLimRob = if_else(Cmpd == 'NSB', Med - 3 * MAD, Med + 3 * MAD)) %>%
  pivot_longer(cols = 6:11, names_to = 'StatType', values_to = 'StatVal') %>%
  pivot_wider(names_from = Cmpd | StatType, names_sep = '_', values_from = StatVal) %>%
  mutate(Z = (NSB_ZLim - TOTB_ZLim) / (NSB_Avg - TOTB_Avg),
         Zrob = (NSB_ZLimRob - TOTB_ZLimRob) / (NSB_Med - TOTB_Med))

ZComp <- PlateQCData %>%
  select(Assay, Scale, AssayPlate, Z, Zrob)%>%
  filter(Scale =='Mean') %>%
  pivot_longer(cols = 4:5, names_to = 'ZFactor', values_to = 'Zval')

ZCompTgt1 <- ZComp %>%
  filter(Assay == 'Tgt1') %>%
  t_test(Zval ~ ZFactor, paired = TRUE) %>%
  add_significance() %>%
  add_xy_position(x = 'ZFactor')

Tgt1Zplot <- ggpaired(filter(ZComp, Assay == 'Tgt1'), x = 'ZFactor', y = 'Zval',
         id = 'AssayPlate',
         order = c('Z', 'Zrob'),
         fill = 'ZFactor',
         ylab = 'Z Value', xlab = 'Z Type') +
  stat_pvalue_manual(ZCompTgt1, tip.length = 0) +
  labs(title = 'Tgt1 Z-Factors',
       subtitle = get_test_label(ZCompTgt1, detailed= TRUE)) +
  theme(legend.position = 'none')

ZCompTgt2 <- ZComp %>%
  filter(Assay == 'Tgt2') %>%
  t_test(Zval ~ ZFactor, paired = TRUE) %>%
  add_significance() %>%
  add_xy_position(x = 'ZFactor')

Tgt2Zplot <- ggpaired(filter(ZComp, Assay == 'Tgt2'), x = 'ZFactor', y = 'Zval',
         id = 'AssayPlate',
         order = c('Z', 'Zrob'),
         fill = 'ZFactor',
         ylab = 'Z Value', xlab = 'Z Type') +
  stat_pvalue_manual(ZCompTgt2, tip.length = 0) +
  labs(title = 'Tgt2 Z-Factors',
       subtitle = get_test_label(ZCompTgt2, detailed= TRUE)) +
  theme(legend.position = 'none')

Fig4 <- (Tgt1Zplot / Tgt2Zplot) +
  plot_annotation(title = 'Figure 4. Standard and Robust Z Comparisons.', tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave('Figures/Weidner Fig 4.jpg', plot = Fig4, height = 6, width = 6, units = 'in', dpi = 300)

# Cmpd Data - remove all wells without test compounds and put into tidy format ----------

CmpdData <- DoseData %>%
  filter(Conc != 0) %>%
  select(-Data, -Row, -Column) %>%
  group_by(Assay, AssayPlate) %>%
  pivot_longer(cols = starts_with('Pct'), names_to = 'Scale', values_to = 'Activity') %>%
  mutate(Scale = if_else(str_detect(Scale, 'Mean'), 'Mean', 'Median')) %>%
  ungroup()

PlateSmplData <- CmpdData %>%
  group_by(Assay, Run, AssayPlate, Sample, Scale) %>%
  mutate(Activity = if_else(Scale == 'Mean', mean(Activity), median(Activity))) %>%
  ungroup()

SampleSummData <- PlateSmplData %>%
  group_by(Assay, Sample, Scale) %>%
  mutate(Activity = if_else(Scale == 'Mean', mean(Activity), median(Activity))) %>%
  ungroup()

# Active Well determination ------------------------
# Working with CmpdData to define "True" activity based on overall mean or median activity per "compound"
# Compound ID is the value in the Sample column
# Calculate overall mean and median for each sample in each assay
# First, split PctActivity into two columns, one for values based on mean controls, one based on median controls

CmpdData2 = CmpdData %>% pivot_wider (names_from = "Scale", values_from = "Activity", names_prefix = "PctAct.")

EstTruth = CmpdData2 %>% group_by (Assay, Sample) %>%
  summarize (true.n = sum (!is.na (PctAct.Mean)),
             true.mean.est = mean (PctAct.Mean),
             true.median.est = median (PctAct.Median))

# Calculate compound well summaries per Assay and plate, using both mean and median, regardless
# of how the control wells are summarized.  In other words, we're going to look at 4 possiblities:
# 1. control wells summarized by mean, compound wells summarized by mean
# 2. control wells summarized by mean, compound wells summarized by median
# 3. control wells summarized by median, compound wells summarized by mean
# 4. control wells summarized by median, compound wells summarized by median

SummPerPlate = CmpdData2 %>% group_by (Assay, AssayPlate, PlateId, Sample) %>%
  summarize (sample.n = sum (!is.na (PctAct.Mean)),
             mean.PctAct.Mean = mean (PctAct.Mean),
             median.PctAct.Mean = median (PctAct.Mean),
             mean.PctAct.Median = mean (PctAct.Median),
             median.PctAct.Median = median (PctAct.Median))

# Merge plate summary data with estimated true values

SummPerPlate2 = merge (SummPerPlate, EstTruth, by=c("Assay", "Sample"))

### Plot the individual well values vs. estimate true values

CmpdData3 = merge (CmpdData2, EstTruth, by=c("Assay", "Sample"))

### Plot summary value per plate vs. estimated true value, means or medians. Fig 5

T1EstTrueMeanN1 <- ggplot (filter(CmpdData3, Assay == 'Tgt1'), aes(x=true.mean.est, y=PctAct.Mean)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  labs(title = 'Tgt1 (N = 1)' )

T2EstTrueMeanN1 <- ggplot (filter(CmpdData3, Assay == 'Tgt2'), aes(x=true.mean.est, y=PctAct.Mean)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  labs(title = 'Tgt2 (N = 1)' )

T1EstTrueMedN1 <- ggplot (filter(CmpdData3, Assay == 'Tgt1'), aes(x=true.mean.est, y=PctAct.Median)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Median)') +
  theme_classic() +
  labs(title = 'Tgt1 (N = 1)' )

T2EstTrueMedN1 <- ggplot (filter(CmpdData3, Assay == 'Tgt2'), aes(x=true.mean.est, y=PctAct.Median)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Median)') +
  theme_classic() +
  labs(title = 'Tgt2 (N = 1)' )

T1EstTrueMeanN4 <- ggplot (filter(SummPerPlate2, Assay == 'Tgt1'), aes(x=true.mean.est, y=mean.PctAct.Mean)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  labs(title = 'Tgt1 (N = 4)' )

T2EstTrueMeanN4 <- ggplot (filter(SummPerPlate2, Assay == 'Tgt2'), aes(x=true.mean.est, y=mean.PctAct.Mean)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  labs(title = 'Tgt2 (N = 4)' )

T1EstTrueMedN4 <- ggplot (filter(SummPerPlate2, Assay == 'Tgt1'), aes(x=true.mean.est, y=median.PctAct.Median)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Median)') +
  theme_classic() +
  labs(title = 'Tgt1 (N = 4)' )

T2EstTrueMedN4 <- ggplot (filter(SummPerPlate2, Assay == 'Tgt2'), aes(x=true.mean.est, y=median.PctAct.Median)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Median)') +
  theme_classic() +
  labs(title = 'Tgt2 (N = 4)' )

T1Scatter <- (T1EstTrueMeanN1/T1EstTrueMedN1/T1EstTrueMeanN4/T1EstTrueMedN4)

T2Scatter <- (T2EstTrueMeanN1/T2EstTrueMedN1/T2EstTrueMeanN4/T2EstTrueMedN4)

Fig5 <- (T1Scatter | T2Scatter) +
  plot_annotation(title = 'Figure 5. Measured Data vs. Estimated Truth', tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave('Figures/Weidner Fig 5.jpg', plot = Fig5, height = 8, width = 6, units = 'in', dpi = 300)

### Function to create ROC curve ------------------------------------------------------
# Uses package, ROCR
# Note that each vertex on the curve corresponds to a particular activity threshold
# Precicted = observed percent activity per compound per plate
# Actual = overall mean or median percent activity per compound
# true.cutoff = cut-off for active using the "actual" values

ROCcurve = function (predicted, actual, true.cutoff = 50) {
  true.labels = ifelse (actual >= true.cutoff, 1, 0)
  pred1 <- prediction(predicted, true.labels)
  perf1 <- performance(pred1,"tpr","fpr")
  auc1 <- performance(pred1,"auc")@y.values[[1]]
  auc1
  plot(perf1, lwd=2, col=2, main = paste ("Activity Cutoff ", true.cutoff, "%", sep=""))
  abline(0,1)
  legend(x = "bottomright", c(paste ("AUC=", round (auc1, 4), sep="")),   lwd=2, col=2)

  # Extract the X and Y values from the ROC plot, as well as the cutoffs
  roc.x = slot (perf1, "x.values") [[1]]
  roc.y = slot (perf1, "y.values") [[1]]
  cutoffs = slot (perf1, "alpha.values") [[1]]

  auc.table = cbind.data.frame(cutoff=pred1@cutoffs,
                               tp=pred1@tp, fp=pred1@fp, tn=pred1@tn, fn=pred1@fn)
  names (auc.table) = c("Cutoff", "TP", "FP", "TN", "FN")
  auc.table$true.cutoff = true.cutoff
  auc.table$sensitivity = auc.table$TP / (auc.table$TP + auc.table$FN)
  auc.table$specificity = auc.table$TN / (auc.table$TN + auc.table$FP)
  auc.table$FalsePosRate = 1 - auc.table$specificity
  auc.table$FalseNegRate = 1 - auc.table$sensitivity
  auc.table$sens_spec = auc.table$sensitivity + auc.table$specificity
  auc.table$PPV = auc.table$TP / (auc.table$TP + auc.table$FP)
  auc.table$NPV = auc.table$TN / (auc.table$TN + auc.table$FN)

  auc.best = auc.table [auc.table$sens_spec == max (auc.table$sens_spec),]
  #row.names (auc.best) = "NULL"

  return (list (roc.table = auc.table, roc.best = auc.best))
}

### Assay Tgt1, means

ROC.mean.mean.50 = with (SummPerPlate2 [SummPerPlate2$Assay == "Tgt1", ],
                         ROCcurve (mean.PctAct.Mean, true.mean.est, 50))
#ROC.mean.mean.50$roc.best

ROC.mean.mean.30 = with (SummPerPlate2 [SummPerPlate2$Assay == "Tgt1", ],
                         ROCcurve (mean.PctAct.Mean, true.mean.est, 30))
ROC.mean.mean.40 = with (SummPerPlate2 [SummPerPlate2$Assay == "Tgt1", ],
                         ROCcurve (mean.PctAct.Mean, true.mean.est, 40))
ROC.mean.mean.60 = with (SummPerPlate2 [SummPerPlate2$Assay == "Tgt1", ],
                         ROCcurve (mean.PctAct.Mean, true.mean.est, 60))
ROC.mean.mean.70 = with (SummPerPlate2 [SummPerPlate2$Assay == "Tgt1", ],
                         ROCcurve (mean.PctAct.Mean, true.mean.est, 70))

ROC.mean.results = rbind (ROC.mean.mean.30$roc.best,
                          ROC.mean.mean.40$roc.best,
                          ROC.mean.mean.50$roc.best,
                          ROC.mean.mean.60$roc.best,
                          ROC.mean.mean.70$roc.best)

## The approach above provides an ROC curve and it chooses an optimal activity cutoff that maximizes
## sensitivity + specificity.  First, a cutoff for true activity is set (30, 40, 50, 60, or 70).
## Then for each of those cutoffs, an optimal cutoff for the plate mean or median is determined.
## I'm not sure this is the best approach, so below I am switching to using the same cutoff for the
## estimated true values and the assay screening values (plate mean or individual wells).  This seems
## to have better balance between false positives and false negatives.  With the ROC approach, there's
## more weight given to false negatives since there is more negative (inactive) data.

### Use the same cutoff for plate results as for the estimated true values

same.cutoff = function (predicted, actual, cutoff) {
  TP = sum (predicted >= cutoff & actual >= cutoff)
  TN = sum (predicted <  cutoff & actual <  cutoff)
  FP = sum (predicted >= cutoff & actual <  cutoff)
  FN = sum (predicted <  cutoff & actual >= cutoff)
  sensitivity = TP / (TP + FN)
  specificity = TN / (TN + FP)
  FalsePosRate = 1 - specificity
  FalseNegRate = 1 - sensitivity
  sens_spec = sensitivity + specificity
  PPV = TP / (TP + FP)
  NPV = TN / (TN + FN)
  return (cbind (cutoff, TP, FP, TN, FN, sensitivity, specificity, FalsePosRate, FalseNegRate, sens_spec, PPV, NPV))
}

PerPlate.Tgt1 = SummPerPlate2 [SummPerPlate2$Assay == "Tgt1", ]

mean.results.Tgt1 = rbind.data.frame (
  with (PerPlate.Tgt1, same.cutoff (mean.PctAct.Mean, true.mean.est, 45)),
  with (PerPlate.Tgt1, same.cutoff (mean.PctAct.Mean, true.mean.est, 50)),
  with (PerPlate.Tgt1, same.cutoff (mean.PctAct.Mean, true.mean.est, 55)),
  with (PerPlate.Tgt1, same.cutoff (mean.PctAct.Mean, true.mean.est, 60)),
  with (PerPlate.Tgt1, same.cutoff (mean.PctAct.Mean, true.mean.est, 65)),
  with (PerPlate.Tgt1, same.cutoff (mean.PctAct.Mean, true.mean.est, 70)),
  with (PerPlate.Tgt1, same.cutoff (mean.PctAct.Mean, true.mean.est, 75)))
mean.results.Tgt1$Scale = "Mean"

### Assay Tgt1, medians

median.results.Tgt1 = rbind.data.frame (
  with (PerPlate.Tgt1, same.cutoff (median.PctAct.Median, true.median.est, 45)),
  with (PerPlate.Tgt1, same.cutoff (median.PctAct.Median, true.median.est, 50)),
  with (PerPlate.Tgt1, same.cutoff (median.PctAct.Median, true.median.est, 55)),
  with (PerPlate.Tgt1, same.cutoff (median.PctAct.Median, true.median.est, 60)),
  with (PerPlate.Tgt1, same.cutoff (median.PctAct.Median, true.median.est, 65)),
  with (PerPlate.Tgt1, same.cutoff (median.PctAct.Median, true.median.est, 70)),
  with (PerPlate.Tgt1, same.cutoff (median.PctAct.Median, true.median.est, 75)))
median.results.Tgt1$Scale = "Median"

### Plot PPV and NPV results, Mean vs Median, for Tgt1

plot.results.Tgt1 = pivot_longer (rbind (mean.results.Tgt1, median.results.Tgt1),
                                  c("PPV", "NPV"), names_to = "Result", values_to = "Value")

# Re-order PPV and NPV

plot.results.Tgt1$Result = factor (plot.results.Tgt1$Result, levels = c("PPV", "NPV"))

T1N4ppv <- ggplot(plot.results.Tgt1, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt1, N=4/Sample") + xlab("Activity Cutoff, %") + ylab("PPV or NPV") +
  theme_classic() +
  theme(legend.position = 'none')

### Plot sensitivity and specificity results, Mean vs Median, for Tgt1

plot.results.Tgt1 = pivot_longer (rbind (mean.results.Tgt1, median.results.Tgt1),
                                  c("sensitivity", "specificity"), names_to = "Result", values_to = "Value")

T1N4sel <- ggplot(plot.results.Tgt1, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt1, N=4/Sample") + xlab("Activity Cutoff, %") + ylab("Sensitivity or Specificity") +
  theme_classic() +
  theme(legend.position = 'none')

### Assay Tgt2, means

PerPlate.Tgt2 = SummPerPlate2 [SummPerPlate2$Assay == "Tgt2", ]

mean.results.Tgt2 = rbind.data.frame (
  with (PerPlate.Tgt2, same.cutoff (mean.PctAct.Mean, true.mean.est, 45)),
  with (PerPlate.Tgt2, same.cutoff (mean.PctAct.Mean, true.mean.est, 50)),
  with (PerPlate.Tgt2, same.cutoff (mean.PctAct.Mean, true.mean.est, 55)),
  with (PerPlate.Tgt2, same.cutoff (mean.PctAct.Mean, true.mean.est, 60)),
  with (PerPlate.Tgt2, same.cutoff (mean.PctAct.Mean, true.mean.est, 65)),
  with (PerPlate.Tgt2, same.cutoff (mean.PctAct.Mean, true.mean.est, 70)),
  with (PerPlate.Tgt2, same.cutoff (mean.PctAct.Mean, true.mean.est, 75)))
mean.results.Tgt2$Scale = "Mean"

### Assay Tgt2, medians

median.results.Tgt2 = rbind.data.frame (
  with (PerPlate.Tgt2, same.cutoff (median.PctAct.Median, true.median.est, 45)),
  with (PerPlate.Tgt2, same.cutoff (median.PctAct.Median, true.median.est, 50)),
  with (PerPlate.Tgt2, same.cutoff (median.PctAct.Median, true.median.est, 55)),
  with (PerPlate.Tgt2, same.cutoff (median.PctAct.Median, true.median.est, 60)),
  with (PerPlate.Tgt2, same.cutoff (median.PctAct.Median, true.median.est, 65)),
  with (PerPlate.Tgt2, same.cutoff (median.PctAct.Median, true.median.est, 70)),
  with (PerPlate.Tgt2, same.cutoff (median.PctAct.Median, true.median.est, 75)))
median.results.Tgt2$Scale = "Median"

### Plot Mean vs Median results for Tgt2

plot.results.Tgt2 = pivot_longer (rbind (mean.results.Tgt2, median.results.Tgt2),
                                  c("PPV", "NPV"), names_to = "Result", values_to = "Value")

# Re-order PPV and NPV

plot.results.Tgt2$Result = factor (plot.results.Tgt2$Result, levels = c("PPV", "NPV"))

T2N4ppv <- ggplot(plot.results.Tgt2, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Target 2, N=4/Sample") + xlab("Activity Cutoff, %") + ylab("PPV or NPV") +
  theme_classic()

### Plot sensitivity and specificity results, Mean vs Median, for Tgt2

plot.results.Tgt2 = pivot_longer (rbind (mean.results.Tgt1, median.results.Tgt2),
                                  c("sensitivity", "specificity"), names_to = "Result", values_to = "Value")

T2N4sel <- ggplot(plot.results.Tgt2, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt2, N=4/Sample") + xlab("Activity Cutoff, %") + ylab("Sensitivity or Specificity") +
  theme_classic()

######################################################################
### Repeat the same-cutoff analysis above for the individual well values

IndivWells.Tgt1 = CmpdData3 [CmpdData3$Assay == "Tgt1", ]

### Assay Tgt1, means

Indiv.meanCtrl.Tgt1 = rbind.data.frame (
  with (IndivWells.Tgt1, same.cutoff (PctAct.Mean, true.mean.est, 45)),
  with (IndivWells.Tgt1, same.cutoff (PctAct.Mean, true.mean.est, 50)),
  with (IndivWells.Tgt1, same.cutoff (PctAct.Mean, true.mean.est, 55)),
  with (IndivWells.Tgt1, same.cutoff (PctAct.Mean, true.mean.est, 60)),
  with (IndivWells.Tgt1, same.cutoff (PctAct.Mean, true.mean.est, 65)),
  with (IndivWells.Tgt1, same.cutoff (PctAct.Mean, true.mean.est, 70)),
  with (IndivWells.Tgt1, same.cutoff (PctAct.Mean, true.mean.est, 75)))
Indiv.meanCtrl.Tgt1$Scale = "Mean"

### Assay Tgt1, medians

Indiv.medianCtrl.Tgt1 = rbind.data.frame (
  with (IndivWells.Tgt1, same.cutoff (PctAct.Median, true.median.est, 45)),
  with (IndivWells.Tgt1, same.cutoff (PctAct.Median, true.median.est, 50)),
  with (IndivWells.Tgt1, same.cutoff (PctAct.Median, true.median.est, 55)),
  with (IndivWells.Tgt1, same.cutoff (PctAct.Median, true.median.est, 60)),
  with (IndivWells.Tgt1, same.cutoff (PctAct.Median, true.median.est, 65)),
  with (IndivWells.Tgt1, same.cutoff (PctAct.Median, true.median.est, 70)),
  with (IndivWells.Tgt1, same.cutoff (PctAct.Median, true.median.est, 75)))
Indiv.medianCtrl.Tgt1$Scale = "Median"

### Plot PPV and NPV results, Mean vs Median, for Tgt1

plot.Indiv.Tgt1 = pivot_longer (rbind (Indiv.meanCtrl.Tgt1, Indiv.medianCtrl.Tgt1),
                                c("PPV", "NPV"), names_to = "Result", values_to = "Value")

# Re-order PPV and NPV

plot.Indiv.Tgt1$Result = factor (plot.Indiv.Tgt1$Result, levels = c("PPV", "NPV"))

T1N1ppv = ggplot(plot.Indiv.Tgt1, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt1, N=1/Sample") + xlab("Activity Cutoff, %") + ylab("PPV or NPV") +
  theme_classic() +
  theme(legend.position = 'none')

### Plot sensitivity and specificity results, Mean vs Median, for Tgt1

plot.Indiv.Tgt1 = pivot_longer (rbind (Indiv.meanCtrl.Tgt1, Indiv.medianCtrl.Tgt1),
                                c("sensitivity", "specificity"), names_to = "Result", values_to = "Value")

T1N1sel <- ggplot(plot.Indiv.Tgt1, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt1, N=1/Sample") + xlab("Activity Cutoff, %") + ylab("Sensitivity or Specificity") +
  theme_classic() +
  theme(legend.position = 'none')

### Assay Tgt2, means

IndivWells.Tgt2 = CmpdData3 [CmpdData3$Assay == "Tgt2", ]

Indiv.meanCtrl.Tgt2 = rbind.data.frame (
  with (IndivWells.Tgt2, same.cutoff (PctAct.Mean, true.mean.est, 45)),
  with (IndivWells.Tgt2, same.cutoff (PctAct.Mean, true.mean.est, 50)),
  with (IndivWells.Tgt2, same.cutoff (PctAct.Mean, true.mean.est, 55)),
  with (IndivWells.Tgt2, same.cutoff (PctAct.Mean, true.mean.est, 60)),
  with (IndivWells.Tgt2, same.cutoff (PctAct.Mean, true.mean.est, 65)),
  with (IndivWells.Tgt2, same.cutoff (PctAct.Mean, true.mean.est, 70)),
  with (IndivWells.Tgt2, same.cutoff (PctAct.Mean, true.mean.est, 75)))
Indiv.meanCtrl.Tgt2$Scale = "Mean"

### Assay Tgt2, medians

Indiv.medianCtrl.Tgt2 = rbind.data.frame (
  with (IndivWells.Tgt2, same.cutoff (PctAct.Median, true.median.est, 45)),
  with (IndivWells.Tgt2, same.cutoff (PctAct.Median, true.median.est, 50)),
  with (IndivWells.Tgt2, same.cutoff (PctAct.Median, true.median.est, 55)),
  with (IndivWells.Tgt2, same.cutoff (PctAct.Median, true.median.est, 60)),
  with (IndivWells.Tgt2, same.cutoff (PctAct.Median, true.median.est, 65)),
  with (IndivWells.Tgt2, same.cutoff (PctAct.Median, true.median.est, 70)),
  with (IndivWells.Tgt2, same.cutoff (PctAct.Median, true.median.est, 75)))
Indiv.medianCtrl.Tgt2$Scale = "Median"

### Plot Mean vs Median results for Tgt1

plot.Indiv.Tgt2 = pivot_longer (rbind (Indiv.meanCtrl.Tgt2, Indiv.medianCtrl.Tgt2),
                                c("PPV", "NPV"), names_to = "Result", values_to = "Value")

# Re-order PPV and NPV

plot.Indiv.Tgt2$Result = factor (plot.Indiv.Tgt2$Result, levels = c("PPV", "NPV"))

T2N1ppv = ggplot(plot.Indiv.Tgt2, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt2, N=1/Sample") + xlab("Activity Cutoff, %") + ylab("PPV or NPV") +
  theme_classic()

### Plot sensitivity and specificity results, Mean vs Median, for Tgt2

plot.Indiv.Tgt2 = pivot_longer (rbind (Indiv.meanCtrl.Tgt1, Indiv.medianCtrl.Tgt2),
                                c("sensitivity", "specificity"), names_to = "Result", values_to = "Value")

T2N1sel <- ggplot(plot.Indiv.Tgt2, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt2, N=1/Sample") + xlab("Activity Cutoff, %") + ylab("Sensitivity or Specificity") +
  theme_classic()

## Put the four plots above together into one plot. Fig 6

T1ROC <- (T1N1sel/T1N1ppv/T1N4sel/T1N4ppv)

T2ROC <- (T2N1sel/T2N1ppv/T2N4sel/T2N4ppv)

Fig6 <- (T1ROC | T2ROC) +
  plot_annotation(title = 'Figure 6. Assay Performance Measures', tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave('Figures/Weidner Fig 6.jpg', plot = Fig6, height = 8, width = 6, units = 'in', dpi = 300)

