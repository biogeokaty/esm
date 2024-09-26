# ESM Data Wrangling for Methods Comparsion

# Setup ----
# Packages
library(here)
library(janitor)
library(readxl)
library(writexl)
library(viridis)
library(paletteer)
library(rms)
library(aqp)
library(cowplot)
library(lme4)
library(tidyverse)

# Data
esm_gen <- read.csv(here("data_processed", "esm_genetic_hrz.csv"))
esm_standard <- read.csv(here("data_processed", "esm_standard_depths.csv"))

# Make data compatible for joining - ESM standard depth increments data ----
# esm_standard has wrong info coded into the "depth increments" column - fix
# also going to remove the repeated fd data
esm_standard <- esm_standard %>%
  mutate(depth_increments = "standard") %>%
  mutate(ref_stat =  case_when(method== "fd" ~ "fd",
                               method!="fd" ~ ref_stat)) %>%
  unite("method_long", c("method", "ref_stat"), sep="_", remove = FALSE)

esm_std_dupes <- esm_standard %>%
  group_by_all() %>%
  filter(n() > 1) %>%
  add_count() # the fixed depth data is duplicated 3 times each (one for each ref stat)

# remove duplicated fd values from esm standard depth increment data
esm_standard_clean <- esm_standard %>%
  group_by_all() %>%
  distinct()

# Make data compatible for joining - ESM genetic horizon depth increments data ----
# remove duplicated fd values from esm genetic horizon depth increment data, add longer method col
esm_gen_clean <- esm_gen %>%
  group_by_all() %>%
  distinct() %>%
  unite("method_long", c("method", "ref_stat"), sep="_", remove = FALSE)

# Join standard and genetic horizon data together
esm_join <- esm_gen_clean %>% 
  rbind(esm_standard_clean)
# Save CSV
write_csv(esm_join, here("data_processed", "esm_join.csv"))