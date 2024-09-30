# Calculating reference soil masses for standard depth increments with reference masses calculated separately for each DSP4SH project

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

# Clean data ----
# Filter out UTRGV and Texas A&M Pt2 data, and any pedons with missing bulk density or SOC 
soc_horizon_filt1 <- soc_horizon %>%
  filter(project!="UTRGV", project!="TexasA&MPt-2") %>% # filter out the projects that are missing data
  group_by(dsp_pedon_id) %>%
  filter(!any(is.na(soc_fill))) %>% # filter out projects missing SOC
  filter(!any(is.na(bd_fill)))  # filter out projects missing BD

# Filter out pedons with max depth < 100 cm 
shallow <- soc_horizon_filt1 %>%
  group_by(dsp_pedon_id) %>%
  filter(hrzdep_b == max(hrzdep_b)) %>%
  filter(hrzdep_b < 100) %>%
  distinct(dsp_pedon_id) %>%
  pull()

soc_horizon_filt <- soc_horizon_filt1 %>%
  filter(!dsp_pedon_id %in% shallow) %>%
  ungroup()

# Generate input spreadsheets for SimpleESM ----

# First, nest the soil horizon data by project
nested <- soc_horizon_filt %>%
  group_by(project) %>%
  mutate(project_copy = project) %>%
  nest()

# Pull out list of projects
project_list <- nested %>%
  pull(project)

# Make nested dataframe with the desired depths
depth_df <- data.frame(
  project = project_list,
  top = rep(c(0,10,30,50,75), each=8),
  bottom = rep(c(10,30,50,75,100), each=8)
) %>%
  arrange(project) %>%
  group_by(project) %>%
  nest() %>%
  rename(depths = data)

# Join soil horizon data and depth dataframe
join <- nest_join(nested, depth_df, by="project") %>%
  unnest(cols=c(depth_df)) %>%
  mutate(depth_list = purrr::map(depths, ~pull(.x, bottom)))

# Run mass aggregation function to get masses for each project
ref_mass_df <- join %>%
  mutate(ref_mass = purrr::map2(.x = data, .y = depth_list, .f = soil_mass_aggregate)) %>%
  unnest(cols=c(depths, ref_mass)) %>%
  select(project, data, top, bottom, mass_agg_min, mass_agg_max, mass_agg_mean) %>%
  mutate(Upper_cm = -top,
         Lower_cm = -bottom) %>%
  select(-top, -bottom) %>%
  nest(ref_mass = c(Upper_cm, Lower_cm, mass_agg_min, mass_agg_max, mass_agg_mean))

# Make the input sheets - grouped by project and summary stat for reference soil mass - each group will have a dataframe for SOC concentrations, BD, and reference soil mass
esm_input <- ref_mass_df %>%
  select(project, data, ref_mass) %>%
  mutate(ref_mass_longer = purrr::map(ref_mass, ~pivot_longer(.x,
                                                              cols=c(mass_agg_min:mass_agg_mean),
                                                              names_to="stat",
                                                              values_to="mass"))) %>%
  # Make the first sheet
  mutate(Concentrations = purrr::map(data, .f = ~{
    select(.x, project_copy, label, dsp_plot_id, dsp_pedon_id, hrzdep_t, hrzdep_b, soc_fill) %>%
      rename("Campaign" = "project_copy",
             "Treatment" = "label",
             "Plot" = "dsp_plot_id",
             "Point" = "dsp_pedon_id",
             "Upper_cm" = "hrzdep_t",
             "Lower_cm" = "hrzdep_b",
             "SOC_g_kg" = "soc_fill") %>%
      mutate(Upper_cm = -Upper_cm,
             Lower_cm = -Lower_cm,
             SOC_g_kg = SOC_g_kg * 10)
  }
  )) %>%
  # Make the second sheet - bulk density
  mutate(BD = purrr::map(data, .f = ~{
    select(.x, project_copy, label, dsp_plot_id, dsp_pedon_id, hrzdep_t, hrzdep_b, bd_fill) %>%
      rename("Campaign" = "project_copy",
             "Treatment" = "label",
             "Plot" = "dsp_plot_id",
             "Point" = "dsp_pedon_id",
             "Upper_cm" = "hrzdep_t",
             "Lower_cm" = "hrzdep_b",
             "BD_g_cm3" = "bd_fill") %>%
      mutate(Upper_cm = -Upper_cm,
             Lower_cm = -Lower_cm)
  })) %>%
  unnest(cols=c(ref_mass_longer)) %>%
  select(project, Upper_cm, Lower_cm, stat, mass, Concentrations, BD) %>%
  group_by(project, stat) %>%
  nest(Ref_soil_mass = c(Upper_cm, Lower_cm, mass)) %>%
  # Make the third sheet - Reference soil mass
  mutate(Ref_soil_mass = purrr::map(Ref_soil_mass, .f = ~{
    select(.x, Upper_cm, Lower_cm, mass) %>%
      mutate(Layer = row_number()) %>%
      rename("Ref_soil_mass_t_ha" = "mass") %>%
      relocate(Layer, Upper_cm, Lower_cm, Ref_soil_mass_t_ha)
  }))

# Write the Excel sheets from the nested dataframe
# make vectors of project/summary stat
project_stat <- esm_input %>% distinct(project, stat)
project_list <- project_stat %>% pull(project) %>% as.character
stat_list <- project_stat %>% pull(stat) %>% as.character

dir.create(here("data_processed", "esm_input_standard"))

map2(.x = project_list,
     .y = stat_list,
     .f = ~{
       conc <- esm_input %>% 
         filter(project==.x, stat==.y) %>%
         ungroup() %>%
         select(Concentrations) %>%
         unnest(cols=c(Concentrations))
       
       bd <- esm_input %>% 
         filter(project==.x, stat==.y) %>%
         ungroup() %>%
         select(BD) %>%
         unnest(cols=c(BD))
       
       mass <- esm_input %>% 
         filter(project==.x, stat==.y) %>%
         ungroup() %>%
         select(Ref_soil_mass) %>%
         unnest(cols=c(Ref_soil_mass))
       
       write_xlsx(list("Concentrations" = conc,
                       "BD" = bd,
                       "Ref_soil_mass" = mass),
                  here("data_processed", "esm_input_standard", glue::glue(.x, "_", .y, ".xlsx")))
     })

# Run SimpleESM function iteratively for all input spreadsheets ----
# Set universal options (these will not change with each input, so we can assign them outside of the map() function)

# Option for the reference soil mass ("manual" or "auto")
RefM_option <- "manual"

# Option for calculations - elements: one or two elements ("SOC_only" or "SOC_and_N")
E_calc_option <- "SOC_only"

# Option for calculations - isotopes: 13C, or 13C and 15N, or not ("13C" or "13C_15N" or "no")
I_calc_option <- "no"

# Create the big folder for the output
dir.create(here("data_processed", "esm_output_standard"))

# Load the SimpleESM function
source(here("code","SimpleESM_function.R"))

# Map function to run the SimpleESM function across all project/summary stat combinations
# This takes a long time to run so be patient! :)
map2(.x = project_list,
     .y = stat_list,
     .f = ~{
       
       # Name of the input .xlsx file
       input_file_name <- here("data_processed", "esm_input_standard", 
                               glue::glue(.x, "_", .y, ".xlsx"))
       
       # Name of the output directory
       output_directory_name <- here("data_processed", "esm_output_standard", glue::glue(.x, "_", .y))
       
       # Create output directory for each project/stat combo
       dir.create(output_directory_name)
       
       SimpleESM(input_file_name, output_directory_name, RefM_option, E_calc_option, I_calc_option)
     })

# Fast way Read in and organize all of the SimpleESM output ----

# Read in all the files - read ESM1 and ESM2 output into different objects because their columns are different
# Ignore the FD output because it calculates based on actual fixed depth. We will re-calculate fixed depth SOC stocks for our desired depth increments manually

# Make a list of the ESM1 output files
esm1_files <- list.files(here("data_processed", "esm_output_standard"), pattern = "ESM\\.csv$", recursive=TRUE, full.names=TRUE)
esm1_files_short <- list.files(here("data_processed", "esm_output_standard"), pattern = "ESM\\.csv$", recursive=TRUE)

# Make a list of the ESM2 output files
esm2_files <- list.files(here("data_processed", "esm_output_standard"), pattern = "ESM2\\.csv$", recursive=TRUE, full.names=TRUE)
esm2_files_short <- list.files(here("data_processed", "esm_output_standard"), pattern = "ESM2\\.csv$", recursive=TRUE)

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
soc_agg_df <- join %>%
  mutate(soc_agg = purrr::map2(data, depth_list, soc_stock_fd)) %>%
  select(project, soc_agg) %>%
  unnest(cols=c(soc_agg)) %>%
  arrange(project, dsp_pedon_id)

# Join dataframes together
esm_join <- esm1_df %>%
  left_join(select(esm2_df, ref_stat, sample_id, topdepth_esm2, depth_esm2, soc_esm2), by=c("sample_id", "ref_stat")) %>%
  left_join(soc_agg_df, by=c("sample_id", "dsp_pedon_id", "layer", "project"))

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
