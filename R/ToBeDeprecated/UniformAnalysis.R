# Uniformity Plate Analysis
# All wells contain the total binding (TOTB) control which generates the maximum raw data signal, but represents vehicle (0% biological activity).

#Set up -------------------------------------------------
library(tidyverse)
library(lubridate)
library(viridis)
library(patchwork)

#rescale_mean -------------------------------------------
rescale_mean <- function(x) {
  x / mean(x)
}

# rescale_med ------------------------------------------------
rescale_med <- function(x) {
  x / median(x)
}

# outlier_robust Calculate outliers using box plot criteria----
outlier_robust <- function(x) {
  x > quantile(x, 0.75) + (1.5 * IQR(x)) | x < quantile(x, 0.25) - (1.5 * IQR(x))
}

# Prepare data -------------------

  UniformData <- read_csv('Data/AllData.csv')%>%
  filter(PlateMap == 'Uniform', PlateType != 'C') %>%
  select(-PlateMap) %>%
  mutate(Assay = factor(Assay, levels = c("Tgt1", "Tgt2")),
         PlateType = as.character(PlateType),
         Well = as_factor(Well),
         PlateId = as_factor(PlateId)) %>%
  group_by(PlateId,ExpTime, Assay, PlateType) %>%
  nest() %>%
  mutate(ExpTime = mdy_hms(ExpTime),
         ExpDate = as_date(ExpTime)) %>%
  unnest(cols= c(data)) %>%
  ungroup() %>%
  group_by(PlateId) %>%
  mutate(MedData = rescale_med(Data)) %>%
  ungroup() %>%
  pivot_longer(cols = contains('Data'),
               names_to = 'Scale',
               values_to = 'Data') %>%
  mutate(Scale = if_else(Scale == 'Data', 'Raw', 'Median'),
         Scale = as_factor(Scale),
         Scale = fct_rev(Scale)
  ) %>%
  ungroup()

T1Raw<- ggplot(filter(UniformData, Assay == 'Tgt1', Scale == 'Raw'), aes(x = PlateId, y = Data, color = PlateType)) +
  geom_boxplot() +
  labs(x = 'Plate',
       y = 'Raw Data') +
  scale_color_viridis_d() +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = 'none') +
  labs(title = 'Tgt1 Raw Data')

T2Raw<- ggplot(filter(UniformData, Assay == 'Tgt2', Scale == 'Raw'), aes(x = PlateId, y = Data, color = PlateType)) +
  geom_boxplot() +
  labs(x = 'Plate',
       y = 'Data') +
  scale_color_viridis_d() +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = 'none') +
  labs(title = 'Tgt2 Raw Data')

T1Med<- ggplot(filter(UniformData, Assay == 'Tgt1', Scale == 'Median'), aes(x = PlateId, y = Data, color = PlateType)) +
  geom_boxplot() +
  labs(x = 'Plate',
       y = 'Median Scaled Data') +
  scale_color_viridis_d() +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = 'none') +
  labs(title = 'Tgt1 Median Scaled')

T2Med<- ggplot(filter(UniformData, Assay == 'Tgt2', Scale == 'Median'), aes(x = PlateId, y = Data, color = PlateType)) +
  geom_boxplot() +
  labs(x = 'Plate',
       y = 'Median Scaled Data') +
  scale_color_viridis_d() +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = 'Tgt2 Median Scaled')

Fig2 <- (T1Raw + T2Raw + T1Med + T2Med) +
  plot_annotation(title = 'Figure 2.', tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave('Figures/Weidner Fig 2.jpg', plot = Fig2, height = 4, width = 8, units = 'in', dpi = 300)
