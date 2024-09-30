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
esm_same <- read.csv(here("data_processed", "esm_same_ref_mass.csv"))

# Join data together ----
esm_gen_std <- esm_gen %>%
  rbind(esm_standard) %>%
  mutate(ref_data = "indv_project")

esm_same2 <- esm_same %>%
  mutate(ref_data = "all_data")

esm_all <- esm_gen_std %>%
  rbind(esm_same2) %>%
  unite("method_long", c("method", "ref_stat"), sep="_", remove = FALSE)

# Save CSV
write_csv(esm_all, here("data_processed", "esm_all.csv"))
