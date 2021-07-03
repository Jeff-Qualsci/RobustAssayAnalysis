# Import plate data files

library(tidyverse)
library(lubridate)

rescale_med <- function(x) {
  x / median(x)
}

PlateCols <- c("PlateOrder","Repeat", "Well", "Type","Data")

# Plate File Parser ----------

InputFiles <- dir('Plates/',  pattern = '*.Txt')
TempData <- tibble() # Initialize an empty dataframe
PlateId <- data.frame(PlateId = 1:length(InputFiles) )

for (i in 1:length(InputFiles)) {
  TimeStamp <- read_tsv(file = paste0('Plates/', InputFiles[i]), skip = 1548, col_names = c('ExpTime'), n_max = 1)

  InFile <- read_tsv(file = paste0('Plates/', InputFiles[i]), skip = 5, col_names = PlateCols, n_max = 1536) %>%
    mutate(File = str_sub(InputFiles[i], 1, -5),
           ExpTime = str_sub(TimeStamp$ExpTime, 45, -1)
    )

  TempData <- rbind(TempData, InFile)
}

TempData <- TempData %>%
  select(File, ExpTime, Well, Data) %>%
  group_by(File, ExpTime) %>%
  mutate(MedScaleData = rescale_med(Data),
         Active = MedScaleData < 0.2,
         PlateMap = if_else(sum(Active) > 100, 'Dose', 'Uniform')) %>%
  select(-MedScaleData, -Active) %>%
  ungroup() %>%
  group_by(File, ExpTime, PlateMap) %>%
  nest() %>%
  mutate(PlateType = if_else(str_detect(File, '3727'), 'A', if_else(str_detect(File, '782'), 'B', 'C')),
         Assay = if_else(str_detect(File, 'dhe') | str_detect(ExpTime, '1/14') | str_detect(ExpTime, '12/10'), 'Tgt1', 'Tgt2'),
         ExpDateTime = mdy_hms(ExpTime)) %>%
  arrange(ExpDateTime) %>%
  bind_cols(PlateId) %>%
  unnest(cols = c(data)) %>%
  ungroup() %>%
  select(PlateId, Assay, ExpTime, PlateMap, PlateType, Well, Data)

write_csv(TempData, 'Data/AllData.csv')
