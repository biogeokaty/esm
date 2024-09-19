# Calculating reference soil masses for standard depth increments

# Setup ----

# Packages
library(here)
library(janitor)
library(readxl)
library(viridis)
library(rms)
library(aqp)
library(tidyverse)

# Data
esm1 <- read_excel(here("data_processed", "Output_ESM.xlsx"), sheet="Output_ESM")
esm2 <- read_excel(here("data_processed", "Output_ESM2.xlsx"), sheet="Output_ESM2")
fd <- read_excel(here("data_processed", "Output_FD.xlsx"), sheet="Output_FD")
soc_pedon <- read.csv(here("data_processed", "04_soc_stock_pedon.csv"))
soc_horizon <- read.csv(here("data_processed", "04_soc_stock_horizon.csv"))

# Functions
# Function to calculate aggregated masses in each depth increment and return max, min, and average values
# Required input: dataframe of soil horizon data, vector with desired bottom depths of each depth increment
soil_mass_aggregate <- function(input, depth){
  # Promote to SPC and dice into 1-cm increments
  spc <- input
  aqp::depths(spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
  
  # dice into 1-cm intervals
  dice <- aqp::dice(spc, fm=0:99 ~ bd_fill)
  
  dice_mass <- horizons(dice) %>%
    mutate(mass = (hrzdep_b - hrzdep_t) * bd_fill * 100)
  
  #Initialize output vector
  out <- vector("list", length(depth))
  
  min_max_mean <- list(
    min = ~min(.x, na.rm=TRUE), 
    max = ~max(.x, na.rm=TRUE),
    mean = ~mean(.x, na.rm=TRUE)
  )
  
  for (i in seq_along(depth)) {
    if (i == 1) {
      out[[i]] <- dice_mass %>%
        filter(hrzdep_t < depth[[i]]) %>%
        group_by(dsp_pedon_id) %>%
        dplyr::summarize(mass_agg = sum(mass)) %>%
        mutate(depth_cat = depth[[i]])
      
    } else {
      out[[i]] <- dice_mass %>%
        filter(hrzdep_t >=depth[[i-1]] & hrzdep_t < depth[[i]]) %>%
        group_by(dsp_pedon_id) %>%
        dplyr::summarize(mass_agg = sum(mass)) %>%
        mutate(depth_cat = depth[[i]])
    }
  }
  
  out_bind <- dplyr::bind_rows(out)
  
  out_bind %>%
    group_by(depth_cat) %>%
    dplyr::summarize(across(mass_agg, min_max_mean))
  
}

# Function to calculate SOC stocks with depth increments that mirror ESM increments
# Required input: dataframe of soil horizon data that includes columns bd_fill, soc_fill, and coarse_frag_fill, vectorwith desired bottom depths of each depth increment
soc_stock_fd <- function(input, depth){
  # Promote to SPC and dice into 1-cm increments
  spc <- input
  aqp::depths(spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
  
  # dice into 1-cm intervals
  dice <- aqp::dice(spc, fm=0:99 ~ bd_fill + soc_fill + coarse_frag_fill)
  
  soc <- horizons(dice) %>%
    mutate(hrzdepth = hrzdep_b - hrzdep_t,
           cf_mult = 1 - (coarse_frag_fill/100)) %>%
    mutate(soc_stock_hrz = soc_fill * bd_fill * hrzdepth * cf_mult)
  
  out <- vector("list", length(depth))
  
  for (i in seq_along(depth)) {
    if (i == 1) {
      out[[i]] <- soc %>%
        filter(hrzdep_t < depth[[i]]) %>%
        group_by(dsp_pedon_id) %>%
        dplyr::summarize(soc_fd = sum(soc_stock_hrz)) %>%
        mutate(depth_cat = depth[[i]])
      
    } else {
      out[[i]] <- soc %>%
        filter(hrzdep_t >=depth[[i-1]] & hrzdep_t < depth[[i]]) %>%
        group_by(dsp_pedon_id) %>%
        dplyr::summarize(soc_fd = sum(soc_stock_hrz)) %>%
        mutate(depth_cat = depth[[i]])
    }
  }
  
  out_bind <- dplyr::bind_rows(out) %>%
    group_by(dsp_pedon_id) %>%
    mutate(layer=seq_along(depth_cat)) %>%
    mutate(depth_fd = -depth_cat) %>%
    select(-depth_cat) %>%
    unite("sample_id", c("dsp_pedon_id", "layer"), sep="-", remove=FALSE) %>%
    ungroup()
  
  out_bind
  
}


# Input data for mass aggregation function ----
# Filter to get the data for the project that you want
project_data <- soc_horizon %>%
  filter(project=="project_name_here")

# Make vector of desired horizon depths
depth_list <- c(10,30,50,75,100)

# Calculated aggregated soil masses and save output ----
mass_agg <- soil_mass_aggregate(project_data, depth_list)
## Use these data to fill in your Reference soil mass sheet for the SimpleESM input spreadsheet! ##

# Calculate fixed depth SOC stocks with increments that mirror ESM depth increments ---
soc_agg <- soc_stock_fd(project_data, depth_list)

# Join ESM calculations and fixed depth calculations into one dataframe ----
# Read in ESM data, add a "sample_id" column for easy joining
esm1 <- read.csv("path_to_ESM1_output_file_goes_here.csv", sep = ";") %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE)
esm2 <- read.csv("path_to_ESM2_output_file_goes_here.csv", sep = ";") %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE)

# Join into one dataframe
soc_join <- esm1 %>%
  select(Campaign, Treatment, Point, sample_id, Layer, Lower_cm, SOC_stock_ESM) %>%
  rename(project = Campaign,
         label=Treatment,
         dsp_pedon_id = Point,
         depth_esm1 = Lower_cm,
         soc_esm1 = SOC_stock_ESM) %>%
  left_join(select(esm2,sample_id, Lower_cm, SOC_stock_ESM2), by="sample_id") %>%
  rename(soc_esm2 = SOC_stock_ESM2,
         depth_esm2 = Lower_cm) %>%
  left_join(soc_agg, by=c("sample_id", "dsp_pedon_id")) %>%
  select(-layer) %>%
  pivot_longer(cols = depth_esm1:depth_fd,
               names_to = c(".value", "method"),
               names_sep="_"
  ) %>%
  rename(soc_stock_calc = soc)

# Generating input spreadsheets for SimpleESM automatically ----
# Requires some fun with nested dataframes
# First, nest the soil horizon data by project
nested <- soc_horizon %>%
  group_by(project) %>%
  mutate(project_copy = project) %>%
  nest()

# Pull out list of projects
project_list <- nested %>%
  pull(project)

# Make nested dataframe with the desired depths
depth_df <- data.frame(
  project = project_list,
  top = rep(c(0,10,30,50,75), each=10),
  bottom = rep(c(10,30,50,75,100), each=10)
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
  mutate(ref_mass = purrr::map2(data, depth_list, soil_mass_aggregate)) %>% 
  unnest(cols=c(depths, ref_mass)) %>%
  select(project, data, top, bottom, mass_agg_min, mass_agg_max, mass_agg_mean) %>%
  mutate(Upper_cm = -top,
         Lower_cm = -bottom) %>%
  select(-top, -bottom) %>%
  nest(ref_mass = c(Upper_cm, Lower_cm, mass_agg_min, mass_agg_max, mass_agg_mean))

# Make the input sheets - grouped by project and summary stat for reference soil mass - each group will have a dataframe for SOC concentrations, BD, and reference soil mass
esm_input <- ref_mass_df %>%
  select(project, data, ref_mass) %>%
  # Filter out any pedons with NA values for SOC fill or BD fill (these will trip up the ESM script)
  mutate(data = purrr::map(data, .f = ~{
    .x %>%
      group_by(dsp_pedon_id) %>%
      filter(!any(is.na(soc_fill))) %>%
      filter(!any(is.na(bd_fill)))
  })) %>%
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
                  here("data_processed", "esm_input_test", glue::glue(.x, "_", .y, ".xlsx")))
     })

# Run SimpleESM function iteratively for all input spreadsheets ----
# Universal options (these will not change with each input, so we can assign them outside of the map() function)

# Option for the reference soil mass ("manual" or "auto")
RefM_option <- "manual"

# Option for calculations - elements: one or two elements ("SOC_only" or "SOC_and_N")
E_calc_option <- "SOC_only"

# Option for calculations - isotopes: 13C, or 13C and 15N, or not ("13C" or "13C_15N" or "no")
I_calc_option <- "no"

# Create the big folder for the output
dir.create(here("data_processed", "esm_output"))

# Load the SimpleESM function
source(here("code","SimpleESM_function.R"))

# Map function to run the SimpleESM function across all project/summary stat combinations
# This takes a long time to run so be patient! :)
map2(.x = project_list,
     .y = stat_list,
     .f = ~{
       
       # Name of the input .xlsx file
       input_file_name <- here("data_processed", "esm_input_test", 
                               glue::glue(.x, "_", .y, ".xlsx"))
       
       # Name of the output directory
       output_directory_name <- here("data_processed", "esm_output", glue::glue(.x, "_", .y))
       
       # Create output directory for each project/stat combo
       dir.create(output_directory_name)
       
       SimpleESM(input_file_name, output_directory_name, RefM_option, E_calc_option, I_calc_option)
     })
