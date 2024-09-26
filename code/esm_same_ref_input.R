# Calculating reference soil masses for standard depth increments with one set of reference masses for all DSP4SH projects

# Setup ----

# Packages
library(here)
library(janitor)
library(readxl)
library(writexl)
library(viridis)
library(rms)
library(aqp)
library(tidyverse)

# Data
soc_horizon <- read.csv(here("data_processed", "04_soc_stock_horizon.csv"))

# Functions
source(here("code","esm_functions.R"))

# Fast way - generating input spreadsheets for SimpleESM automatically ----
# Calculate min, mean, and max soil mass in 0-10 cm, 10-30 cm, 30-50 cm, 50-75 cm, and 75-10 cm depths 
depths <- c(10,30,50,75,100)

soc_horizon_filt <- soc_horizon %>%
  filter(project!="UTRGV", project!="TexasA&MPt-2") # filter out the projects that are missing data
 
ref_mass <- soil_mass_aggregate(soc_horizon_filt, depths)

# Build SOC and BD input dataframes ----
esm_input_soc <- soc_horizon_filt %>%
  group_by(dsp_pedon_id) %>%
  filter(!any(is.na(soc_fill))) %>%
  filter(!any(is.na(bd_fill))) %>%
  select(project, label, dsp_plot_id, dsp_pedon_id, hrzdep_t, hrzdep_b, soc_fill) %>%
  rename("Campaign" = "project",
         "Treatment" = "label",
         "Plot" = "dsp_plot_id",
         "Point" = "dsp_pedon_id",
         "Upper_cm" = "hrzdep_t",
         "Lower_cm" = "hrzdep_b",
         "SOC_g_kg" = "soc_fill") %>%
  mutate(Upper_cm = -Upper_cm,
         Lower_cm = -Lower_cm,
         SOC_g_kg = SOC_g_kg * 10)

esm_input_bd <- soc_horizon_filt %>%
  group_by(dsp_pedon_id) %>%
  filter(!any(is.na(soc_fill))) %>%
  filter(!any(is.na(bd_fill))) %>%
  select(project, label, dsp_plot_id, dsp_pedon_id, hrzdep_t, hrzdep_b, bd_fill) %>%
  rename("Campaign" = "project",
         "Treatment" = "label",
         "Plot" = "dsp_plot_id",
         "Point" = "dsp_pedon_id",
         "Upper_cm" = "hrzdep_t",
         "Lower_cm" = "hrzdep_b",
         "BD_g_cm3" = "bd_fill") %>%
  mutate(Upper_cm = -Upper_cm,
         Lower_cm = -Lower_cm)
  
  
# Make dataframe with the desired depths
depth_df <- data.frame(
  top = c(0,10,30,50,75),
  bottom = c(10,30,50,75,100)
)

# Join depth dataframe with aggregated soil masses
mass_agg_df <- ref_mass %>%
  left_join(depth_df, by=c("depth_cat" = "bottom")) %>%
  rename("bottom" = "depth_cat") %>%
  mutate(Layer = row_number(),
         top = -top,
         bottom = -bottom) %>%
  pivot_longer(mass_agg_min:mass_agg_mean, 
         names_to="stat",
         values_to="mass") %>%
  group_by(stat) %>%
  nest() %>%
  mutate(data = purrr::map(data, .f = ~{
    select(.x, Layer, bottom, top, mass) %>%
      rename("Ref_soil_mass_t_ha" = "mass",
             "Upper_cm" = "top",
             "Lower_cm" = "bottom") %>%
      relocate(Layer, Upper_cm, Lower_cm, Ref_soil_mass_t_ha)
  }))
  
# Write the Excel spreadsheets
stat_list <- mass_agg_df %>% pull(stat) %>% as.character

purrr::map(.x = stat_list,
     .f = ~{
       
       esm_input_mass <- mass_agg_df %>% 
         filter(stat== .x) %>%
         ungroup() %>%
         select(data) %>%
         unnest(cols=c(data))
       
       write_xlsx(list("Concentrations" = esm_input_soc,
                       "BD" = esm_input_bd,
                       "Ref_soil_mass" = esm_input_mass),
                  here("data_processed", "esm_input_same_ref_mass", glue::glue("esm_input_same_ref", .x, ".xlsx")))
     })

# Fast way - Run SimpleESM function iteratively for all input spreadsheets ----
# Set universal options (these will not change with each input, so we can assign them outside of the map() function)

# Option for the reference soil mass ("manual" or "auto")
RefM_option <- "manual"

# Option for calculations - elements: one or two elements ("SOC_only" or "SOC_and_N")
E_calc_option <- "SOC_only"

# Option for calculations - isotopes: 13C, or 13C and 15N, or not ("13C" or "13C_15N" or "no")
I_calc_option <- "no"

# Create the big folder for the output
dir.create(here("data_processed", "esm_output_same_ref_mass"))

# Load the SimpleESM function
source(here("code","SimpleESM_function.R"))

# Map function to run the SimpleESM function across all project/summary stat combinations
# This takes a long time to run so be patient! :)
map(.x = stat_list,
     .f = ~{
       # Name of the input .xlsx file
       input_file_name <- here("data_processed", "esm_input_same_ref_mass", 
                               glue::glue("esm_input_same_ref", .x, ".xlsx"))
       
       # Name of the output directory
       output_directory_name <- here("data_processed", "esm_output_same_ref_mass", glue::glue(.x))
       
       # Create output directory for each project/stat combo
       dir.create(output_directory_name)
       
       SimpleESM(input_file_name, output_directory_name, RefM_option, E_calc_option, I_calc_option)
     })

# Read in and organize all of the SimpleESM output ----

# Read in all the files - read ESM1 and ESM2 output into different objects because their columns are different
# Ignore the FD output because it calculates based on actual fixed depth. We will re-calculate fixed depth SOC stocks for our desired depth increments manually

# Make a list of the ESM1 output files
esm1_files <- list.files(here("data_processed", "esm_output_same_ref_mass"), pattern = "ESM\\.csv$", recursive=TRUE, full.names=TRUE)
esm1_files_short <- list.files(here("data_processed", "esm_output_same_ref_mass"), pattern = "ESM\\.csv$", recursive=TRUE)

# Make a list of the ESM2 output files
esm2_files <- list.files(here("data_processed", "esm_output_same_ref_mass"), pattern = "ESM2\\.csv$", recursive=TRUE, full.names=TRUE)
esm2_files_short <- list.files(here("data_processed", "esm_output_same_ref_mass"), pattern = "ESM2\\.csv$", recursive=TRUE)

# Read in files and name each one so we can tell them apart
esm1_data <- lapply(esm1_files, read.csv, sep = ";")
names(esm1_data) <- gsub("\\.csv$", "", esm1_files_short)

esm2_data <- lapply(esm2_files, read.csv, sep = ";")
names(esm2_data) <- gsub("\\.csv$", "", esm2_files_short)

# Convert the lists of dataframes into actual dataframes
# Start with ESM1 dataframe
esm1_df <- do.call(rbind.data.frame, esm1_data) %>%
  rownames_to_column() %>%
  mutate(ref_stat = case_when(grepl("max", rowname) ~ "max",
                              grepl("min", rowname) ~ "min",
                              grepl("mean", rowname) ~ "mean")) %>%
  select(-rowname) %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE) %>%
  select(Campaign, ref_stat, Treatment, Point, sample_id, Layer, Upper_cm, Lower_cm, SOC_stock_ESM) %>%
  rename(project = Campaign,
         label=Treatment,
         dsp_pedon_id = Point,
         topdepth_esm1 = Upper_cm,
         depth_esm1 = Lower_cm,
         soc_esm1 = SOC_stock_ESM,
         layer = Layer)

# Convert ESM2 list of data into dataframe
esm2_df <- do.call(rbind.data.frame, esm2_data) %>%
  rownames_to_column() %>%
  mutate(ref_stat = case_when(grepl("max", rowname) ~ "max",
                              grepl("min", rowname) ~ "min",
                              grepl("mean", rowname) ~ "mean")) %>%
  select(-rowname) %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE) %>%
  select(Campaign, ref_stat, Treatment, Point, sample_id, Layer, Upper_cm, Lower_cm, SOC_stock_ESM2) %>%
  rename(project = Campaign,
         label=Treatment,
         dsp_pedon_id = Point,
         topdepth_esm2 = Upper_cm,
         depth_esm2 = Lower_cm,
         soc_esm2 = SOC_stock_ESM2,
         layer = Layer)

# Calculate fixed depth SOC stocks for all projects to join in
soc_agg_df <- soc_stock_fd(soc_horizon_filt, depths)

# Join dataframes together
esm_join <- esm1_df %>%
  left_join(select(esm2_df, ref_stat, sample_id, topdepth_esm2, depth_esm2, soc_esm2), by=c("sample_id", "ref_stat")) %>%
  left_join(soc_agg_df, by=c("sample_id", "dsp_pedon_id", "layer"))

# Make longer
esm_join_long <- esm_join %>%
  pivot_longer(cols = topdepth_esm1:depth_fd,
               names_to = c(".value", "method"),
               names_sep="_") %>%
  mutate(depth_increments = "standard") %>%
  mutate(apparent_depth = case_when(method == "fd" & layer==1 ~ glue::glue("0{round(depth, 1)} cm"),
                                    method == "fd" & layer > 1 ~ glue::glue("{-round(topdepth, 1)}{round(depth,1)} cm")),
         actual_depth = case_when(layer==1 ~ glue::glue("0{round(depth, 1)} cm"),
                                  layer > 1 ~ glue::glue("{-round(topdepth, 1)}{round(depth,1)} cm"))) %>%
  group_by(project, layer) %>%
  fill(apparent_depth, .direction="up") %>%
  mutate(ref_stat =  case_when(method== "fd" ~ "fd",
                               method!="fd" ~ ref_stat)) %>%
  group_by_all() %>%
  distinct()

# Write csv
write_csv(esm_join_long, here("data_processed", "esm_standard_depths.csv"))
