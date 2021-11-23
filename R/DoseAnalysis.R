# Dose Response Data Analysis
# All wells contain the total binding (TOTB) control which generates the maximum raw data signal, but represents vehicle (0% biological activity).

#Set up -------------------------------------------------
library(tidyverse)
library(lubridate)
library(viridis)
library(ggpubr)
library(rstatix)
#library(ROCR)
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

# Create Control Well Plots Fig 3------
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

# Comparison of Standard and Robust Z-factors ------------
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
# Define "True" activity based on overall mean or median activity per sample

CmpdData <- DoseData %>%
  filter(Conc != 0) %>%
  select(-Data, -Row, -Column) %>%
  group_by(Assay, AssayPlate) %>%
  pivot_longer(cols = starts_with('Pct'), names_to = 'Scale', values_to = 'Activity') %>%
  mutate(Scale = if_else(str_detect(Scale, 'Mean'), 'Mean', 'Median')) %>%
  ungroup() %>%
  group_by(Assay, Sample, Scale) %>%
  mutate(SampTrueActivity = if_else(Scale == 'Mean', mean(Activity), median(Activity))) %>%
  ungroup()

# PlateSmplData - Summarize 4 plate replicates and then determine estimated truth ----------
# Plate summarization will be consistent with scaling method Median%Act summarized as medians
# Mean%Act summarized as means
# Active Well determination,define "True" activity based on overall mean or median activity per sample, 4 samples/plate

PlateSmplData <- CmpdData %>%
  group_by(Assay, Run, AssayPlate, Sample, Scale) %>%
  mutate(PlateSampleActivity = if_else(Scale == 'Mean', mean(Activity), median(Activity))) %>%
  select(-Well, -Activity, -SampTrueActivity) %>%
  ungroup() %>%
  unique() %>%
  group_by(Assay, Sample, Scale) %>%
  mutate(PlateSampTrueActivity = if_else(Scale == 'Mean', mean(PlateSampleActivity), median(PlateSampleActivity))) %>%
  ungroup()

### Plot summary value per plate vs. estimated true value, means or medians. Fig 5 -----------

T1EstTrueMeanN1 <- ggplot (filter(CmpdData, Assay == 'Tgt1', Scale == 'Mean'), aes(x=SampTrueActivity, y=Activity)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  labs(title = 'Tgt1 (N = 1)' )

T2EstTrueMeanN1 <- ggplot (filter(CmpdData, Assay == 'Tgt2', Scale == 'Mean'), aes(x=SampTrueActivity, y=Activity)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  labs(title = 'Tgt2 (N = 1)' )

T1EstTrueMedN1 <- ggplot (filter(CmpdData, Assay == 'Tgt1', Scale == 'Median'), aes(x=SampTrueActivity, y=Activity)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Median)') +
  theme_classic() +
  labs(title = 'Tgt1 (N = 1)' )

T2EstTrueMedN1 <- ggplot (filter(CmpdData, Assay == 'Tgt2', Scale == 'Median'), aes(x=SampTrueActivity, y=Activity)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Median)') +
  theme_classic() +
  labs(title = 'Tgt2 (N = 1)' )

T1EstTrueMeanN4 <- ggplot (filter(PlateSmplData, Assay == 'Tgt1', Scale == 'Mean'), aes(x = PlateSampTrueActivity, y = PlateSampleActivity)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  labs(title = 'Tgt1 (N = 4)' )

T2EstTrueMeanN4 <- ggplot (filter(PlateSmplData, Assay == 'Tgt2', Scale == 'Mean'), aes(x = PlateSampTrueActivity, y = PlateSampleActivity)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  labs(title = 'Tgt2 (N = 4)' )

T1EstTrueMedN4 <- ggplot (filter(PlateSmplData, Assay == 'Tgt1', Scale == 'Median'), aes(x = PlateSampTrueActivity, y = PlateSampleActivity)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Median)') +
  theme_classic() +
  labs(title = 'Tgt1 (N = 4)' )

T2EstTrueMedN4 <- ggplot (filter(PlateSmplData, Assay == 'Tgt2', Scale == 'Median'), aes(x = PlateSampTrueActivity, y = PlateSampleActivity)) +
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

# Assign Estimated and measured Active values based on various threshold values. ----------

CmpdData <- CmpdData %>%
  mutate(TrueAct_30 = SampTrueActivity >= 30,
         MeasAct_30 = Activity >= 30,
         TrueAct_40 = SampTrueActivity >= 40,
         MeasAct_40 = Activity >= 40,
         TrueAct_50 = SampTrueActivity >= 50,
         MeasAct_50 = Activity >= 60,
         TrueAct_60 = SampTrueActivity >= 60,
         MeasAct_60 = Activity >= 60,
         TrueAct_70 = SampTrueActivity >= 70,
         MeasAct_70 = Activity >= 70,
         TrueAct_80 = SampTrueActivity >= 80,
         MeasAct_80 = Activity >= 80)

# Prepare Sample data for ROC analysis -----------
SampleActData  <-  CmpdData %>%
  select(-PlateType, -ExpTime, -ExpDate, -Run) %>%
  pivot_longer(cols = ends_with('0'),
               names_to = c('ActType', 'Threshold'),
               names_sep = '_',
               names_transform = list(Threshold = as.integer),
               values_to = 'Active')%>%
  pivot_wider(names_from = ActType,
              values_from = Active) %>%
  mutate(TrueAct = unlist(TrueAct),
         MeasAct = unlist(MeasAct),
         ConfusionType = if_else(TrueAct,
                                 if_else(MeasAct,'TP', 'FN'),
                                 if_else(MeasAct, 'FP', 'TN'))
         )

SampleConfusionData <- SampleActData %>%
  group_by(Assay, Scale, Threshold, ConfusionType) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = ConfusionType,
              values_from = Count) %>%
  mutate(Sensitivity = TP / (TP + FN),
         Specificity = TN / (TN + FP),
         PPV = TP / (TP + FP),
         NPV = TN / (TN + FN))

# Prepare Plate Sample data for ROC Analysis -------------
PlateSmplData <- PlateSmplData %>%
  mutate(TrueAct_30 = PlateSampTrueActivity >= 30,
         MeasAct_30 = PlateSampleActivity >= 30,
         TrueAct_40 = PlateSampTrueActivity >= 40,
         MeasAct_40 = PlateSampleActivity >= 40,
         TrueAct_50 = PlateSampTrueActivity >= 50,
         MeasAct_50 = PlateSampleActivity >= 60,
         TrueAct_60 = PlateSampTrueActivity >= 60,
         MeasAct_60 = PlateSampleActivity >= 60,
         TrueAct_70 = PlateSampTrueActivity >= 70,
         MeasAct_70 = PlateSampleActivity >= 70,
         TrueAct_80 = PlateSampTrueActivity >= 80,
         MeasAct_80 = PlateSampleActivity >= 80)

PlateActData  <-  CmpdData %>%
  select(-PlateType, -ExpTime, -ExpDate, -Run) %>%
  pivot_longer(cols = ends_with('0'),
               names_to = c('ActType', 'Threshold'),
               names_sep = '_',
               names_transform = list(Threshold = as.integer),
               values_to = 'Active')%>%
  pivot_wider(names_from = ActType,
              values_from = Active) %>%
  mutate(TrueAct = unlist(TrueAct),
         MeasAct = unlist(MeasAct),
         ConfusionType = if_else(TrueAct,
                                 if_else(MeasAct,'TP', 'FN'),
                                 if_else(MeasAct, 'FP', 'TN'))
  )

PlateConfusionData <- PlateActData %>%
  group_by(Assay, Scale, Threshold, ConfusionType) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = ConfusionType,
              values_from = Count) %>%
  mutate(Sensitivity = TP / (TP + FN),
         Specificity = TN / (TN + FP),
         PPV = TP / (TP + FP),
         NPV = TN / (TN + FN))


# Creat ROC graphs for Figure 6 ----------------
SampleConfusionPPV <- SampleConfusionData %>%
  select(Assay,Scale, Threshold, ends_with('V')) %>%
  pivot_longer(cols = 4:5, names_to = 'PerfMetric', values_to = 'value')

T1N1ppv = ggplot(filter(SampleConfusionPPV, Assay == 'Tgt1'), aes (x=Threshold, y=value, color=PerfMetric, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt1, N=1/Sample") + xlab("Activity Cutoff, %") + ylab("PPV or NPV") +
  theme_classic() +
  theme(legend.position = 'none')

T2N1ppv = ggplot(filter(SampleConfusionPPV, Assay == 'Tgt2'), aes (x=Threshold, y=value, color=PerfMetric, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt2, N=1/Sample") + xlab("Activity Cutoff, %") + ylab("PPV or NPV") +
  theme_classic()

SampleConfusionSen <- SampleConfusionData %>%
  select(Assay,Scale, Threshold, starts_with('S')) %>%
  pivot_longer(cols = 4:5, names_to = 'PerfMetric', values_to = 'value')

T1N1sel <- ggplot(filter(SampleConfusionSen, Assay == 'Tgt1'), aes (x=Threshold, y=value, color=PerfMetric, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.7, 1) +
  ggtitle ("Tgt1, N=1/Sample") + xlab("Activity Cutoff, %") + ylab("Sensitivity or Specificity") +
  theme_classic() +
  theme(legend.position = 'none')

T2N1sel <- ggplot(filter(SampleConfusionSen, Assay == 'Tgt2'), aes (x=Threshold, y=value, color=PerfMetric, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.7, 1) +
  ggtitle ("Tgt2, N=1/Sample") + xlab("Activity Cutoff, %") + ylab("Sensitivity or Specificity") +
  theme_classic()

PlateConfusionPPV <- PlateConfusionData %>%
  select(Assay,Scale, Threshold, ends_with('V')) %>%
  pivot_longer(cols = 4:5, names_to = 'PerfMetric', values_to = 'value')

T1N4ppv = ggplot(filter(PlateConfusionPPV, Assay == 'Tgt1'), aes (x=Threshold, y=value, color=PerfMetric, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt1, N=4/Sample") + xlab("Activity Cutoff, %") + ylab("PPV or NPV") +
  theme_classic() +
  theme(legend.position = 'none')

T2N4ppv = ggplot(filter(PlateConfusionPPV, Assay == 'Tgt2'), aes (x=Threshold, y=value, color=PerfMetric, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt2, N=4/Sample") + xlab("Activity Cutoff, %") + ylab("PPV or NPV") +
  theme_classic()

PlateConfusionSen <- PlateConfusionData %>%
  select(Assay,Scale, Threshold, starts_with('S')) %>%
  pivot_longer(cols = 4:5, names_to = 'PerfMetric', values_to = 'value')

T1N4sel <- ggplot(filter(SampleConfusionSen, Assay == 'Tgt1'), aes (x=Threshold, y=value, color=PerfMetric, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.7, 1) +
  ggtitle ("Tgt1, N=4/Sample") + xlab("Activity Cutoff, %") + ylab("Sensitivity or Specificity") +
  theme_classic() +
  theme(legend.position = 'none')

T2N4sel <- ggplot(filter(SampleConfusionSen, Assay == 'Tgt2'), aes (x=Threshold, y=value, color=PerfMetric, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.7, 1) +
  ggtitle ("Tgt2, N=4/Sample") + xlab("Activity Cutoff, %") + ylab("Sensitivity or Specificity") +
  theme_classic()

# Put the four plots above together into one plot. Fig 6

T1ROC <- (T1N1sel/T1N1ppv/T1N4sel/T1N4ppv)

T2ROC <- (T2N1sel/T2N1ppv/T2N4sel/T2N4ppv)

Fig6 <- (T1ROC | T2ROC) +
  plot_annotation(title = 'Figure 6. Assay Performance Measures', tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave('Figures/Weidner Fig 6.jpg', plot = Fig6, height = 8, width = 6, units = 'in', dpi = 300)

