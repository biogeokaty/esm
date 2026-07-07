# Calculating reference soil masses for standard depth increments with reference masses calculated separately for each DSP4SH project and hybrid bulk density

# Setup ----

# Packages
library(here)
library(janitor)
library(readxl)
library(writexl)
library(rms)
library(aqp)
library(tidyverse)

# Data
soc_horizon <- read.csv(here("data_processed", "04_soc_stock_horizon.csv"))
project <- read.csv(here("data_processed", "project_data_meta.csv"))

# Functions
source(here("code","esm_functions.R"))

# Clean data for ESM calculations ----
# Filter out Texas A&M Pt2 data, and any pedons with missing bulk density or SOC 
soc_horizon_filt1 <- soc_horizon %>%
  filter(project!="TexasA&MPt-2") %>% # filter out the projects that are missing data (just Texas A&M pt 2)
  group_by(dsp_pedon_id) %>%
  filter(!any(is.na(soc_fill))) %>% # filter out pedons missing SOC
  filter(!any(is.na(bd_fill))) %>%  # filter out pedons missing BD
  mutate(project = case_when(project=="OregonState" & soil=="Jory" ~ "OregonStateJory",
                             project=="OregonState" & soil=="Woodburn" ~ "OregonStateWoodburn",
                             .default = project))

# Filter out pedons with max depth < 75 cm 
shallow <- soc_horizon_filt1 %>%
  group_by(dsp_pedon_id) %>%
  filter(hrzdep_b == max(hrzdep_b)) %>%
  filter(hrzdep_b < 75) %>%
  distinct(dsp_pedon_id) %>%
  pull()

# filter out shallow pedons and correct bulk density for coarse fragment content
soc_horizon_filt <- soc_horizon_filt1 %>%
  filter(!dsp_pedon_id %in% shallow) %>%
  ungroup() %>%
  mutate(bd_hybrid = bd_fill * cf_mult) %>%  
  select(project, label, dsp_plot_id, dsp_pedon_id, dsp_sample_id, layer_no, hzdesg, hrzdep_t, hrzdep_b, soc_fill, bd_hybrid)

# write_csv(soc_horizon_filt, here("data_processed", "soc_horizon_esm_input.csv"))

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
  top = rep(c(0,5,10,30), each=10),
  bottom = rep(c(5,10,30,60), each=10)) %>%
  arrange(project) %>%
  group_by(project) %>%
  nest() %>%
  rename(depths = data)

# Join soil horizon data and depth dataframe
join <- nest_join(nested, depth_df, by="project") %>%
  unnest(cols=c(depth_df)) %>%
  mutate(depth_list = purrr::map(depths, ~pull(.x, bottom)))

## Run mass aggregation function to get masses for each project ----

# One version - reference mass data is from the min and mean mass for the entire project
ref_mass_df <- join %>%
  mutate(ref_mass = purrr::map2(.x = data, .y = depth_list, .f = soil_mass_aggregate_hybrid)) %>%
  unnest(cols=c(depths, ref_mass)) %>%
  select(project, data, top, bottom, mass_agg_min, mass_agg_mean, mass_agg_max) %>%
  mutate(Upper_cm = -top,
         Lower_cm = -bottom) %>%
  select(-top, -bottom) %>%
  nest(ref_mass = c(Upper_cm, Lower_cm, mass_agg_min, mass_agg_mean, mass_agg_max))

# Make another version - reference mass data is from the min and mean mass for the treatment within each project expected to have the lowest BD - Reference
ref_mass_refonly_df <- join %>%
  mutate(data_refonly = purrr::map(data, ~filter(.x, label=="Ref")),
         ref_mass = purrr::map2(.x = data_refonly, .y = depth_list, .f = soil_mass_aggregate_hybrid)) %>%
  unnest(cols=c(depths, ref_mass)) %>%
  select(project, data, top, bottom, mass_agg_min, mass_agg_mean, mass_agg_max) %>%
  mutate(Upper_cm = -top,
         Lower_cm = -bottom) %>%
  select(-top, -bottom) %>%
  nest(ref_mass = c(Upper_cm, Lower_cm, mass_agg_min, mass_agg_mean, mass_agg_max))

## Make the input sheets ---- 
# grouped by project and summary stat for reference soil mass - each group will have a dataframe for SOC concentrations, BD, and reference soil mass
### Reference mass from any treatment in project ----
esm_input <- ref_mass_df %>%
  select(project, data, ref_mass) %>%
  mutate(ref_mass_longer = purrr::map(ref_mass, ~pivot_longer(.x,
                                                              cols=c(mass_agg_min, mass_agg_mean, mass_agg_max),
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
    select(.x, project_copy, label, dsp_plot_id, dsp_pedon_id, hrzdep_t, hrzdep_b, bd_hybrid) %>%
      rename("Campaign" = "project_copy",
             "Treatment" = "label",
             "Plot" = "dsp_plot_id",
             "Point" = "dsp_pedon_id",
             "Upper_cm" = "hrzdep_t",
             "Lower_cm" = "hrzdep_b",
             "BD_g_cm3" = "bd_hybrid") %>%
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

### Reference mass from reference treatment only ----
esm_refonly_input <- ref_mass_refonly_df %>%
  select(project, data, ref_mass) %>%
  mutate(ref_mass_longer = purrr::map(ref_mass, ~pivot_longer(.x,
                                                              cols=c(mass_agg_min, mass_agg_mean, mass_agg_max),
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
    select(.x, project_copy, label, dsp_plot_id, dsp_pedon_id, hrzdep_t, hrzdep_b, bd_hybrid) %>%
      rename("Campaign" = "project_copy",
             "Treatment" = "label",
             "Plot" = "dsp_plot_id",
             "Point" = "dsp_pedon_id",
             "Upper_cm" = "hrzdep_t",
             "Lower_cm" = "hrzdep_b",
             "BD_g_cm3" = "bd_hybrid") %>%
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

## Write the Excel sheets from the nested dataframe ----
# make vectors of project/summary stat
project_stat <- esm_input %>% distinct(project, stat)
project_list <- project_stat %>% pull(project) %>% as.character
stat_list <- project_stat %>% pull(stat) %>% as.character

### Reference mass from any treatment in project ----
dir.create(here("data_processed", "esm_revision3"))
dir.create(here("data_processed", "esm_revision3", "esm_input")) # don't need to run this if folder already exists

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
                  here("data_processed", "esm_revision3", "esm_input", glue::glue(.x, "_", .y, ".xlsx")))
     })

### Reference mass from reference treatment only ----
project_stat2 <- esm_refonly_input %>% distinct(project, stat)
project_list2 <- project_stat2 %>% pull(project) %>% as.character
stat_list2 <- project_stat2 %>% pull(stat) %>% as.character

dir.create(here("data_processed", "esm_revision3", "esm_input_refonly"))

map2(.x = project_list2,
     .y = stat_list2,
     .f = ~{
       conc <- esm_refonly_input %>% 
         filter(project==.x, stat==.y) %>%
         ungroup() %>%
         select(Concentrations) %>%
         unnest(cols=c(Concentrations))
       
       bd <- esm_refonly_input %>% 
         filter(project==.x, stat==.y) %>%
         ungroup() %>%
         select(BD) %>%
         unnest(cols=c(BD))
       
       mass <- esm_refonly_input %>% 
         filter(project==.x, stat==.y) %>%
         ungroup() %>%
         select(Ref_soil_mass) %>%
         unnest(cols=c(Ref_soil_mass))
       
       write_xlsx(list("Concentrations" = conc,
                       "BD" = bd,
                       "Ref_soil_mass" = mass),
                  here("data_processed", "esm_revision3", "esm_input_refonly", glue::glue(.x, "_", .y, ".xlsx")))
     })

# Run SimpleESM function iteratively for all input spreadsheets ----
# Set universal options (these will not change with each input, so we can assign them outside of the map() function)

# Option for the reference soil mass ("manual" or "auto")
RefM_option <- "manual"

# Option for calculations - elements: one or two elements ("SOC_only" or "SOC_and_N")
E_calc_option <- "SOC_only"

# Option for calculations - isotopes: 13C, or 13C and 15N, or not ("13C" or "13C_15N" or "no")
I_calc_option <- "no"

### Reference mass from any treatment in project ----
# Create the big folder for the output
dir.create(here("data_processed", "esm_revision3", "esm_output"))

# Load the SimpleESM function
source(here("code","SimpleESM_function.R"))

# Map function to run the SimpleESM function across all project/summary stat combinations
# This takes a long time to run so be patient! :)
map2(.x = project_list,
     .y = stat_list,
     .f = ~{
       
       # Name of the input .xlsx file
       input_file_name <- here("data_processed", "esm_revision3", "esm_input", 
                               glue::glue(.x, "_", .y, ".xlsx"))
       
       # Name of the output directory
       output_directory_name <- here("data_processed", "esm_revision3","esm_output", glue::glue(.x, "_", .y))
       
       # Create output directory for each project/stat combo
       dir.create(output_directory_name)
       
       SimpleESM(input_file_name, output_directory_name, RefM_option, E_calc_option, I_calc_option)
     })

### Reference mass from reference treatment only ----

# Create the big folder for the output
dir.create(here("data_processed", "esm_revision3","esm_output_refonly"))

map2(.x = project_list2,
     .y = stat_list2,
     .f = ~{
       
       # Name of the input .xlsx file
       input_file_name <- here("data_processed", "esm_revision3", "esm_input_refonly", 
                               glue::glue(.x, "_", .y, ".xlsx"))
       
       # Name of the output directory
       output_directory_name <- here("data_processed", "esm_revision3", "esm_output_refonly", glue::glue(.x, "_", .y))
       
       # Create output directory for each project/stat combo
       dir.create(output_directory_name)
       
       SimpleESM(input_file_name, output_directory_name, RefM_option, E_calc_option, I_calc_option)
     })

# Read in and organize all of the SimpleESM output ----

## Reference mass from any treatment in project ----

# Make a list of the ESM1 output files
esm1_files <- list.files(here("data_processed", "esm_revision3", "esm_output"), pattern = "ESM\\.csv$", recursive=TRUE, full.names=TRUE)
esm1_files_short <- list.files(here("data_processed", "esm_revision3", "esm_output"), pattern = "ESM\\.csv$", recursive=TRUE)

# Make a list of the ESM2 output files
esm1_data <- lapply(esm1_files, read.csv, sep = ";")
names(esm1_data) <- gsub("\\.csv$", "", esm1_files_short)

esm2_files <- list.files(here("data_processed", "esm_revision3", "esm_output"), pattern = "ESM2\\.csv$", recursive=TRUE, full.names=TRUE)
esm2_files_short <- list.files(here("data_processed", "esm_revision3", "esm_output"), pattern = "ESM2\\.csv$", recursive=TRUE)

# Read in files and name each one so we can tell them apart
esm2_data <- lapply(esm2_files, read.csv, sep = ";")
names(esm2_data) <- gsub("\\.csv$", "", esm2_files_short)

# Convert the lists of dataframes into actual dataframes
# Start with ESM1 dataframe
esm1_df <- do.call(rbind.data.frame, esm1_data) %>%
  rownames_to_column() %>%
  mutate(ref_stat = case_when(grepl("min", rowname) ~ "min",
                              grepl("mean", rowname) ~ "mean",
                              grepl("max", rowname) ~ "max"),
         ref_data = "project") %>%
  select(-rowname) %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE) %>%
  select(Campaign, ref_stat, ref_data, Treatment, Point, sample_id, Layer, Upper_cm, Lower_cm, 
         Soil_mass_cum_ESM, SOC_stock_cum_ESM, Soil_mass_ESM, SOC_stock_ESM) %>%
  rename(project = Campaign,
         label=Treatment,
         dsp_pedon_id = Point,
         topdepth.esm1 = Upper_cm,
         depth.esm1 = Lower_cm,
         soil_mass_cum.esm1 = Soil_mass_cum_ESM,
         soc_cum.esm1 = SOC_stock_cum_ESM,
         soil_mass.esm1 = Soil_mass_ESM,
         soc.esm1 = SOC_stock_ESM,
         layer = Layer)

# Convert ESM2 list of data into dataframe
esm2_df <- do.call(rbind.data.frame, esm2_data) %>%
  rownames_to_column() %>%
  mutate(ref_stat = case_when(grepl("min", rowname) ~ "min",
                              grepl("mean", rowname) ~ "mean",
                              grepl("max", rowname) ~ "max"),
         ref_data = "project") %>%
  select(-rowname) %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE) %>%
  select(Campaign, ref_stat, ref_data, Treatment, Point, sample_id, Layer, Upper_cm, Lower_cm, 
         Soil_mass_cum_ESM2, SOC_stock_cum_ESM2, Soil_mass_ESM2, SOC_stock_ESM2) %>%
  rename(project = Campaign,
         label=Treatment,
         dsp_pedon_id = Point,
         topdepth.esm2 = Upper_cm,
         depth.esm2 = Lower_cm,
         soil_mass_cum.esm2 = Soil_mass_cum_ESM2,
         soc_cum.esm2 = SOC_stock_cum_ESM2,
         soil_mass.esm2 = Soil_mass_ESM2,
         soc.esm2 = SOC_stock_ESM2,
         layer = Layer)

## Reference mass from reference treatment only ----
# Make a list of the ESM1 output files
esm1_refonly_files <- list.files(here("data_processed", "esm_revision3", "esm_output_refonly"), 
                                 pattern = "ESM\\.csv$", recursive=TRUE, full.names=TRUE)
esm1_refonly_files_short <- list.files(here("data_processed", "esm_revision3", "esm_output_refonly"), 
                                       pattern = "ESM\\.csv$", recursive=TRUE)

# Make a list of the ESM2 output files
esm2_refonly_files <- list.files(here("data_processed", "esm_revision3", "esm_output_refonly"), 
                                 pattern = "ESM2\\.csv$", recursive=TRUE, full.names=TRUE)
esm2_refonly_files_short <- list.files(here("data_processed", "esm_revision3", "esm_output_refonly"), 
                                       pattern = "ESM2\\.csv$", recursive=TRUE)

# Read in files and name each one so we can tell them apart
esm1_refonly_data <- lapply(esm1_refonly_files, read.csv, sep = ";")
names(esm1_refonly_data) <- gsub("\\.csv$", "", esm1_refonly_files_short)

esm2_refonly_data <- lapply(esm2_refonly_files, read.csv, sep = ";")
names(esm2_refonly_data) <- gsub("\\.csv$", "", esm2_refonly_files_short)

# Convert the lists of dataframes into actual dataframes
# Start with ESM1 dataframe
esm1_refonly_df <- do.call(rbind.data.frame, esm1_refonly_data) %>%
  rownames_to_column() %>%
  mutate(ref_stat = case_when(grepl("min", rowname) ~ "min",
                              grepl("mean", rowname) ~ "mean",
                              grepl("max", rowname) ~ "max"),
         ref_data = "treatment") %>%
  select(-rowname) %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE) %>%
  select(Campaign, ref_stat, ref_data, Treatment, Point, sample_id, Layer, Upper_cm, Lower_cm, 
         Soil_mass_cum_ESM, SOC_stock_cum_ESM, Soil_mass_ESM, SOC_stock_ESM) %>%
  rename(project = Campaign,
         label=Treatment,
         dsp_pedon_id = Point,
         topdepth.esm1 = Upper_cm,
         depth.esm1 = Lower_cm,
         soil_mass_cum.esm1 = Soil_mass_cum_ESM,
         soc_cum.esm1 = SOC_stock_cum_ESM,
         soil_mass.esm1 = Soil_mass_ESM,
         soc.esm1 = SOC_stock_ESM,
         layer = Layer)

# Convert ESM2 list of data into dataframe
esm2_refonly_df <- do.call(rbind.data.frame, esm2_refonly_data) %>%
  rownames_to_column() |> 
  mutate(ref_stat = case_when(grepl("min", rowname) ~ "min",
                              grepl("mean", rowname) ~ "mean",
                              grepl("max", rowname) ~ "max"),
         ref_data = "treatment") %>%
  select(-rowname) %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE) %>%
  select(Campaign, ref_stat, ref_data, Treatment, Point, sample_id, Layer, Upper_cm, Lower_cm, 
         Soil_mass_cum_ESM2, SOC_stock_cum_ESM2, Soil_mass_ESM2, SOC_stock_ESM2) %>%
  rename(project = Campaign,
         label=Treatment,
         dsp_pedon_id = Point,
         topdepth.esm2 = Upper_cm,
         depth.esm2 = Lower_cm,
         soil_mass_cum.esm2 = Soil_mass_cum_ESM2,
         soc_cum.esm2 = SOC_stock_cum_ESM2,
         soil_mass.esm2 = Soil_mass_ESM2,
         soc.esm2 = SOC_stock_ESM2,
         layer = Layer)

## Calculate fixed depth SOC and mass ----
# Calculate fixed depth SOC stocks for all projects to join in
# Using this instead of the SimpleESM output because I want SOC stocks calculated in the same fixed depth increments as the ESM depths
# SimpleESM just outputs FD data based on original input data (which is in genetic horizons)
soc_agg_df <- join %>%
  mutate(soc_agg = purrr::map2(data, depth_list, soc_stock_fd_hybrid)) %>%
  select(project, soc_agg) %>%
  unnest(cols=c(soc_agg)) %>%
  arrange(project, dsp_pedon_id) %>%
  pivot_longer(cols=c(topdepth.fd, depth.fd, soc.fd, soil_mass.fd, soc_cum.fd, soil_mass_cum.fd),
               names_to = c(".value", "method"),
               names_pattern="^(.*)\\.(.*)$") %>%
  left_join(select(project, dsp_pedon_id, label), by="dsp_pedon_id") %>% # join in treatment
  mutate(ref_stat="fd",
         ref_data="fd") %>%
  relocate(project, ref_stat, ref_data, label, dsp_pedon_id, sample_id, layer, method, topdepth, depth, 
           soil_mass_cum, soc_cum, soil_mass, soc)

# Join ESM dataframes together and pivot longer ----
esm_project_join <- esm1_df %>%
  left_join(select(esm2_df, ref_stat, ref_data, sample_id, topdepth.esm2, depth.esm2, 
                   soc.esm2, soil_mass.esm2, soc_cum.esm2, soil_mass_cum.esm2), 
            by=c("sample_id", "ref_stat", "ref_data"))

esm_refonly_join <- esm1_refonly_df %>%
  left_join(select(esm2_refonly_df, ref_stat, ref_data, sample_id, topdepth.esm2, depth.esm2, 
                   soc.esm2, soil_mass.esm2, soc_cum.esm2, soil_mass_cum.esm2), 
            by=c("sample_id", "ref_stat", "ref_data"))

esm_join <- bind_rows(esm_project_join, esm_refonly_join) %>%
  pivot_longer(cols = !project:layer,
               names_to = c(".value", "method"),
               names_pattern="^(.*)\\.(.*)$") %>%
  bind_rows(soc_agg_df) %>%
  mutate(esm_depth = case_when(method == "fd" & layer==1 ~ glue::glue("0{round(depth, 1)} cm"),
                               method == "fd" & layer > 1 ~ glue::glue("{-round(topdepth, 1)}{round(depth,1)} cm")),
         actual_depth = case_when(layer==1 ~ glue::glue("0{round(depth, 1)} cm"),
                                  layer > 1 ~ glue::glue("{-round(topdepth, 1)}{round(depth,1)} cm"))) %>%
  group_by(project, layer) %>%
  fill(esm_depth, .direction="up") %>%
  group_by_all() %>%
  distinct()

# Check that # of observations are consistent between methods
esm_join %>% ungroup() %>% count(project, method, ref_stat, ref_data)

# Write csv
write_csv(esm_join, here("data_processed", "esm_revision3", "esm_standard_depths.csv"))
