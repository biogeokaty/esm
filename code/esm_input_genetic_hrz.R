## Generating input data for SimpleESM function, running function, and organizing output
## Scenario: different reference soils for each project, depth increments calculated by genetic horizon

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

# Calculate horizon depth boundaries and reference soil masses for each horizon in each project ----

# Nest soil profile data by project
nested <- soc_horizon_filt %>%
  group_by(project) %>%
  mutate(project_copy = project) %>%
  nest()

# Manually set reference depth masses - want all projects to have a horizon from 0-10 cm, and then the rest of the profile is based on the results of the regression model (see archived section at end of script for where these depths came from)

kansas_depth_df <- data.frame(project = "KansasState",
                              top = c(0, 10, 23, 56, 85),
                              bottom = c(10, 23, 56, 85, 100))
ncs_depth_df <- data.frame(project = "NCState",
                           top = c(0, 10, 16, 28, 57),
                           bottom = c(10, 16, 28, 57, 100))
wash_depth_df <- data.frame(project = "WashingtonState",
                            top = c(0, 10, 30, 59, 87),
                            bottom = c(10, 30, 59, 87, 100))
minn_depth_df <- data.frame(project = "UnivOfMinnesota",
                            top = c(0, 10, 38, 53),
                            bottom = c(10, 38, 53, 100))
ill_depth_df <- data.frame(project= "Illinois",
                           top = c(0, 10, 35, 69),
                           bottom = c(10, 35, 69, 100))
uconn_depth_df <- data.frame(project="UConn",
                             top = c(0, 10, 23, 59, 82),
                             bottom = c(10, 23, 59, 82, 100))
tam1_depth_df <- data.frame(project= "TexasA&MPt-1",
                            top = c(0, 10, 35, 75),
                            bottom = c(10, 35, 75, 100))
osu_depth_df <- data.frame(project= "OregonState",
                           top = c(0, 10, 29, 60),
                           bottom = c(10, 29, 60, 100))

# Bind depths into a nested dataframe
depths_manual <- rbind(kansas_depth_df, ncs_depth_df, wash_depth_df, minn_depth_df, ill_depth_df, 
                       uconn_depth_df, tam1_depth_df, osu_depth_df) %>%
  group_by(project) %>%
  nest() %>%
  rename(depths = data)

# Join depth data with soil data, add vector of depth list
depths_df <- nest_join(nested, depths_manual, by="project") %>%
  select(project, data, depths_manual) %>%
  unnest(cols=c(depths_manual)) %>%
  mutate(depth_list = purrr::map(depths, ~pull(.x, bottom)))

# Calculate reference soil masses
ref_mass_df <- depths_df %>%
  mutate(ref_mass = purrr::map2(data, depth_list, soil_mass_aggregate)) %>% 
  unnest(cols=c(depths, ref_mass)) %>%
  select(project, data, top, bottom, mass_agg_min, mass_agg_max, mass_agg_mean) %>%
  mutate(Upper_cm = -top,
         Lower_cm = -bottom) %>%
  select(-top, -bottom) %>%
  nest(ref_mass = c(Upper_cm, Lower_cm, mass_agg_min, mass_agg_max, mass_agg_mean))

# Make input spreadsheets for SimpleESM function ----

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

dir.create(here("data_processed", "esm_input_genetic"))

# Only run this line when you need to actually write the Excel sheets
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
                  here("data_processed", "esm_input_genetic", glue::glue(.x, "_", .y, ".xlsx")))
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
dir.create(here("data_processed", "esm_output_genetic"))

# Load the SimpleESM function
source(here("code","SimpleESM_function.R"))

# Map function to run the SimpleESM function
map2(.x = project_list,
     .y = stat_list,
     .f = ~{
       
       # Name of the input .xlsx file
       input_file_name <- here("data_processed", "esm_input_genetic", 
                               glue::glue(.x, "_", .y, ".xlsx"))
       
       # Name of the output directory
       output_directory_name <- here("data_processed", "esm_output_genetic", glue::glue(.x, "_", .y))
       
       # Create output directory for each project/stat combo
       dir.create(output_directory_name)
       
       SimpleESM(input_file_name, output_directory_name, RefM_option, E_calc_option, I_calc_option)
     })

# Read in and organize all of the SimpleESM output ----
# Map function to read in all the files - read ESM1 and ESM2 output into different objects because their columns are different
esm1_files <- list.files(here("data_processed", "esm_output_genetic"), pattern = "ESM\\.csv$", recursive=TRUE, full.names=TRUE)
esm1_files_short <- list.files(here("data_processed", "esm_output_genetic"), pattern = "ESM\\.csv$", recursive=TRUE)

esm2_files <- list.files(here("data_processed", "esm_output_genetic"), pattern = "ESM2\\.csv$", recursive=TRUE, full.names=TRUE)
esm2_files_short <- list.files(here("data_processed", "esm_output_genetic"), pattern = "ESM2\\.csv$", recursive=TRUE)

esm1_data <- lapply(esm1_files, read.csv, sep = ";")
names(esm1_data) <- gsub("\\.csv$", "", esm1_files_short)

esm2_data <- lapply(esm2_files, read.csv, sep = ";")
names(esm2_data) <- gsub("\\.csv$", "", esm2_files_short)

# Convert list of dataframes into actual dataframes
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

# Need to also calculate fixed depth SOC stocks for all projects to join in
soc_agg_df <- depths_df %>%
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
  mutate(depth_increments = "genetic_hrz") %>%
  mutate(apparent_depth = case_when(method == "fd" & layer==1 ~ glue::glue("0{round(depth, 1)} cm"),
                                    method == "fd" & layer > 1 ~ glue::glue("{-round(topdepth, 1)}{round(depth,1)} cm")),
         actual_depth = case_when(layer==1 ~ glue::glue("0{round(depth, 1)} cm"),
                                  layer > 1 ~ glue::glue("{-round(topdepth, 1)}{round(depth,1)} cm"))) %>%
  mutate(ref_stat =  case_when(method== "fd" ~ "fd",
                               method!="fd" ~ ref_stat)) %>%
  group_by(project, layer) %>%
  fill(apparent_depth, .direction="up") %>%
  mutate(ref_stat =  case_when(method== "fd" ~ "fd",
                               method!="fd" ~ ref_stat)) %>%
  group_by_all() %>%
  distinct() # filter out fixed depth duplicates

# Write csv
write_csv(esm_join_long, here("data_processed", "esm_genetic_hrz.csv"))

# ARCHIVE - Automatic calculation of depth intervals ----
# Assign generalized horizon labels for all projects ----

# This needs to be done manually by looking at soil profiles, there really is no shortcut
# Kansas
n_ks <- c('A', 'Bt', 'Btk', 'Bk') # generalized horizon label sequence
p_ks <- c('^A',
          '^Bt$|^Bt1|^Bt2|^Bt3|^Bt4',
          '^Btk',
          '^Bk') # regex pattern to assign generalized labels to specific horizons

# NC State
n_ncs <- c('A', 'B', 'Bt1','Bt2','Bt3', 'BC', 'C')
p_ncs <- c('^Ap|^A$',
           '^B$|BA|B/A|A/B|E',
           '^Bt1|^Bt$',
           '^Bt2',
           '^Bt3',
           '^BC|^Bt/BC',
           '^C|\\dC|3Abp|3Bbt')

# TAM1
n_tam1 <- c('A', 'Bt1', 'Bt2') 
p_tam1 <- c('^Ap|^A$',
            'Bt1',
            'Bt2')

# Washington State
n_wash <- c('Ap', "A",'AB', 'Bw') 
p_wash <- c('^Ap',
            '^A$',
            'AB',
            '^Bw')

# University of Minnesota
n_minn <- c('Ap1', 'Ap2', 'A', 'AB', 'BA', 'Bw', 'Bw2', 'Bt') 
p_minn <- c('Ap1',
            'Ap2',
            '^A$|Ap3',
            'AB$',
            'BA$',
            '^Bw|2C1|2BC',
            '^2Bw|2BG',
            'Bt')

# Oregon State
n_osu <- c('A', 'AB', 'Bt')
p_osu <- c('^Ap|^A$',
           'AB|BA',
           '^Bt|\\dBt|2BC|BCt')

# Illinois 
n_ill <- c('A', 'Bt1', 'Bt2', 'C') 
p_ill <- c('^Ap|^A$|AB',
           'Bt1|BAt|Btg1',
           'Bt2|Btg2|2 Bt|Btg3',
           'C')

# UConn
n_uconn <- c('A', 'Bw', 'BC', 'C') 
p_uconn <- c('A',
             '^Bw|Bw$',
             '^BC',
             '^C')

# Put generalized horizon labels and patterns together into list
n_p <- tibble(project=c("KansasState", "NCState", "WashingtonState", "UnivOfMinnesota", 
                        "Illinois", "UConn", "TexasA&MPt-1", "OregonState"),
              name=list(n_ks, n_ncs, n_wash, n_minn, n_ill, n_uconn, n_tam1, n_osu),
              pattern=list(p_ks, p_ncs, p_wash, p_minn, p_ill, p_uconn, p_tam1, p_osu))

# Calculate depth intervals
# Excluding Illinois and UTRGV projects for now because the profiles aren't working - may have to calculate manually

dd <- NULL # need to run this or depth assignment function will not work

depths_df_fun <- join %>%
  filter(project!="Illinois", project!="UTRGV") %>%
  mutate(depths = purrr::pmap(list(data, name, pattern), genhz_depths)) %>%
  select(project, data, depths) %>%
  mutate(depth_list = purrr::map(depths, ~pull(.x, bottom)))

# Looking at the depths we get - it is probably preferable to calculate all depth intervals manually. Not all of the results make sense (e.g. Oregon State with only one horizon boundary at 29 cm)

# Depths to review manually - Illinois (algorithm didn't work), OSU (worked but not very many horizons designated), and UTRGV (algorithm didn't work)
# Illinois

illinois <- soc_horizon %>%
  filter(project=="Illinois") %>%
  filter(hzdesg!="")
spc_illinois <- illinois
depths(spc_illinois) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b #make data a SoilProfileCollection, get depths from hrzdep_t and hrzdep_b columns
hzdesgnname(spc_illinois) <- 'hzdesg' #get horizon names from hzdesg column
plotSPC(spc_illinois, color="hzdesg")

# Dice the profile first so that we don't have to include unnecessary horizons 
dice_illinois <- aqp::dice(spc_illinois, seq(from = 0, to = 99, by = 1) ~ hzdesg)
# Add generalized horizon labels, remove unused levels, and make factor
dice_illinois$genhz <- aqp::generalize.hz(dice_illinois$hzdesg, n_ill, p_ill) 
dice_illinois$genhz[dice_illinois$genhz == "not-used"] <- NA
dice_illinois$genhz <- factor(dice_illinois$genhz)

plotSPC(spc_illinois, color="hzdesg")
plotSPC(dice_illinois, color="genhz")

# keep track of generalized horizon names for later
hz_names_illinois <- levels(dice_illinois$genhz)

# Logistic Proportional-Odds Ordinal Regression Model
horizons_illinois <- aqp::horizons(dice_illinois)
dd_illinois <- rms::datadist(horizons_illinois)
options(datadist = "dd_illinois")
l.genhz_illinois <- orm(genhz ~ rcs(hrzdep_t), data = horizons_illinois, x = TRUE, y = TRUE)

# predict along same depths: columns are the class-wise probability fitted.ind --> return all probability estimates
predict_illinois <- data.frame(predict(l.genhz_illinois, data.frame(hrzdep_t = seq(from = 0, to = 99, by = 1)), type = "fitted.ind"))

# re-name, rms model output give funky names
names(predict_illinois) <- hz_names_illinois

# add depths
predict_illinois$top <- seq(from = 0, to = 99, by = 1)

# maximum likelihood depths
predict_ml_illinois <- get.ml.hz(predict_illinois, o.names = hz_names_illinois)

# return output
predict_ml_illinois

# Oregon State
osu <- soc_horizon %>%
  filter(project=="OregonState")
spc_osu <- osu
depths(spc_osu) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b #make data a SoilProfileCollection, get depths from hrzdep_t and hrzdep_b columns
hzdesgnname(spc_osu) <- 'hzdesg' #get horizon names from hzdesg column
plotSPC(spc_osu, color="hzdesg")

# Add generalized horizon labels
spc_osu$genhz <- aqp::generalize.hz(spc_osu$hzdesg, n_osu, p_osu) 
plotSPC(spc_osu, color="genhz")
# It looks like the reason that the regression model only assigns two horizons is that profiles are mostly an A horizon with a deep Bt - could either re-assign horizon designations (to Bt2, Bt3, etc),  or just eyeball a place in the profile to split