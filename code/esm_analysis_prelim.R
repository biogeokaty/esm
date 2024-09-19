# Initial exploration of ESM data

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
esm1 <- read_excel(here("data_processed", "Output_ESM.xlsx"), sheet="Output_ESM")
esm2 <- read_excel(here("data_processed", "Output_ESM2.xlsx"), sheet="Output_ESM2")
fd <- read_excel(here("data_processed", "Output_FD.xlsx"), sheet="Output_FD")
soc_pedon <- read.csv(here("data_processed", "04_soc_stock_pedon.csv"))
soc_horizon <- read.csv(here("data_processed", "04_soc_stock_horizon.csv"))

# Functions
source(here("code","esm_functions.R"))

# Initial data checks ----
# This dataset: calculated ESM and FD SOC stocks where entire DSP4SH dataset is used as input and SimpleESM calculates reference masses automatically

# want to check that FD gives same results as our fixed depth calculations

fd_100 <- fd %>%
  filter(Lower_cm == -100) %>%
  select(Point, SOC_stock_cum_FD)

# plot script output vs our output
fd_compare <- left_join(soc_pedon, fd_100, by=c("dsp_pedon_id" = "Point"))
ggplot(fd_compare, aes(soc_stock_100cm, SOC_stock_cum_FD)) +
  geom_point() +
  geom_abline(slope=1, intercept=0)

# filter out the observations that don't match
fd_nomatch <- fd_compare %>%
  mutate(ratio = SOC_stock_cum_FD / soc_stock_100cm) %>%
  filter(ratio>1.05) %>%
  pull(dsp_pedon_id)

fd_nomatch_data <- fd %>%
  filter(Point %in% fd_nomatch)

soc_horizon_nomatch <- soc_horizon %>%
  filter(dsp_pedon_id %in% fd_nomatch)

# the samples that don't match are the ones with coarse frags :)

# compare FD at 30 cm to ESM methods at 30 cm ----
# I am assuming that the ESM calculations go to 30 cm?? it's actually really hard to know...

# Pull out the deepest layer for ESM calculations to get cumulative SOC stock
esm1_deepest <- esm1 %>%
  group_by(Campaign, Treatment, Point) %>%
  filter(Lower_cm == min(Lower_cm)) %>%
  ungroup()

esm2_deepest <- esm2 %>%
  group_by(Campaign, Treatment, Point) %>%
  filter(Lower_cm == min(Lower_cm)) %>%
  ungroup()

# Make table with fixed depth SOC stocks and ESM SOC stocks
soc_stock30_compare <- soc_pedon %>%
  select(project, dsp_pedon_id, label, soc_stock_0_30cm) %>%
  rename(stock_fd = soc_stock_0_30cm) %>%
  mutate(depth_fd = -30)%>% 
  left_join(select(esm1_deepest, Point, Lower_cm, SOC_stock_cum_ESM), by=c("dsp_pedon_id" = "Point")) %>%
  rename(stock_esm1 = SOC_stock_cum_ESM,
         depth_esm1 = Lower_cm) %>%
  left_join(select(esm2_deepest, Point, Lower_cm, SOC_stock_cum_ESM2), by=c("dsp_pedon_id" = "Point")) %>%
  rename(stock_esm2 = SOC_stock_cum_ESM2,
         depth_esm2 = Lower_cm) %>%
  pivot_longer(cols = stock_fd:stock_esm2,
               names_to = c(".value", "method"),
               names_sep="_"
               )
  
# Plot comparison of stocks
ggplot(soc_stock30_compare, aes(x=label, y=stock, fill=method)) +
  geom_boxplot() +
  facet_wrap(~project, scales="free") +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()

# Plot comparison of depths used
ggplot(soc_stock30_compare, aes(x=label, y=depth, fill=method)) +
  geom_boxplot() +
  facet_wrap(~project, scales="free") +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()

max(soc_stock30_compare$depth)
# The shallowest depth used by ESM to calculate total SOC stocks is 10 cm. I think this is an artifact of comparing very different soils to each other, rather than reflecting real differences in bulk density driven by management

# Also, looking at the Texas A&M Pt2 data - they didn't measure SOC below 22 cm in most cases - I wonder if the inclusion of those shallow samples is the reason that the ESM calculations are really shallow.

# Try re-running the calculations without the Texas A&M Pt 2 data ----
#  Read in output of SimpleESM script with no TAM2 data

fd_no_tam2 <- read.csv(here("data_processed", "simple_esm", "Output_FD_no_tam2.csv"), sep = ";")
esm1_no_tam2 <- read.csv(here("data_processed", "simple_esm", "Output_ESM_no_tam2.csv"), sep = ";")
esm2_no_tam2 <- read.csv(here("data_processed", "simple_esm", "Output_ESM2_no_tam2.csv"), sep = ";")

min(esm1_no_tam2$Lower_cm)
# getting a little deeper - this time the minimum depth is 37 cm

# ESM calcs for Oregon State project alone ----
fd_osu <- read.csv(here("data_processed", "simple_esm", "Output_FD_osu.csv"), sep = ";")
esm1_osu <- read.csv(here("data_processed", "simple_esm", "Output_ESM_osu.csv"), sep = ";")
esm2_osu <- read.csv(here("data_processed", "simple_esm", "Output_ESM2_osu.csv"), sep = ";")

min(esm1_osu$Lower_cm)
# Minimum depth is -30.2 - seems like the script defaults to only going down to 30 cm if auto method is used?
# I can't figure out why that would be the case from looking at the script
# Why are there only 3 layers?

# Compare our FD calcs to ESM calcs
esm1_deepest_osu <- esm1_osu %>%
  group_by(Campaign, Treatment, Point) %>%
  filter(Lower_cm == min(Lower_cm)) %>%
  ungroup()

esm2_deepest_osu <- esm2_osu %>%
  group_by(Campaign, Treatment, Point) %>%
  filter(Lower_cm == min(Lower_cm)) %>%
  ungroup()

compare_osu <- soc_pedon %>%
  filter(project=="OregonState") %>%
  select(project, dsp_pedon_id, label, soc_stock_0_30cm) %>%
  rename(stock_fd = soc_stock_0_30cm) %>%
  mutate(depth_fd = -30)%>% 
  left_join(select(esm1_deepest_osu, Point, Lower_cm, SOC_stock_cum_ESM), by=c("dsp_pedon_id" = "Point")) %>%
  rename(stock_esm1 = SOC_stock_cum_ESM,
         depth_esm1 = Lower_cm) %>%
  left_join(select(esm2_deepest_osu, Point, Lower_cm, SOC_stock_cum_ESM2), by=c("dsp_pedon_id" = "Point")) %>%
  rename(stock_esm2 = SOC_stock_cum_ESM2,
         depth_esm2 = Lower_cm) %>%
  pivot_longer(cols = stock_fd:stock_esm2,
               names_to = c(".value", "method"),
               names_sep="_"
  )

# Plot comparison of stocks
ggplot(compare_osu, aes(x=label, y=stock, fill=method)) +
  geom_boxplot() +
  labs(y="SOC stock 0-30 cm (Mg/ha)") +
  facet_wrap(~project, scales="free") +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()

# Plot comparison of depths used
ggplot(compare_osu, aes(x=label, y=depth, fill=method)) +
  geom_boxplot() +
  facet_wrap(~project, scales="free") +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()

# Calculate soil mass in different layers for OSU data ----
# Reference depths to use:
# 0-10 cm
# 10-30 cm
# 30-50 cm
# 50-75 cm
# 75-100 cm


# First need to know what a reasonable soil mass is for the desired depths
# Extract data from the SoilProfileCollection back into a dataframe and calculate soil mass in each 1-cm increment
spc_osu <- soc_horizon %>% filter(project=="OregonState")
depths(spc_osu) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b #make data a SoilProfileCollection, get depths from hrzdep_t and hrzdep_b columns
hzdesgnname(spc_osu) <- 'hzdesg' #get horizon names from hzdesg column

# Use dice() function to resample profile into 1 cm increments, keeping the variables needed for SOC stock calculation
dice_osu <- aqp::dice(spc_osu, fm=0:100 ~ soc_fill + bd_fill + coarse_frag_fill)

# calculate mass of each 1-cm horizon
mass_osu <- horizons(dice_osu) %>%
  mutate(mass = (hrzdep_b - hrzdep_t) * bd_fill * 100) # 1 g/cm2 = 100 megagrams ha

# 0-10 cm depth sum
mass0_10_osu <- mass_osu %>%
  filter(hrzdep_t < 10) %>%
  group_by(dsp_pedon_id) %>%
  summarize(mass_tot10 = sum(mass))

# 10-30 cm depth sum
mass10_30_osu <- mass_osu %>%
  filter(hrzdep_t >=10 & hrzdep_t < 30) %>%
  group_by(dsp_pedon_id) %>%
  summarize(mass_tot30 = sum(mass))

# 30-50 cm depth sum
mass30_50_osu <- mass_osu %>%
  filter(hrzdep_t >=30 & hrzdep_t < 50) %>%
  group_by(dsp_pedon_id) %>%
  summarize(mass_tot50 = sum(mass))

# 50-75 cm depth sum
mass50_75_osu <- mass_osu %>%
  filter(hrzdep_t >=50 & hrzdep_t < 75) %>%
  group_by(dsp_pedon_id) %>%
  summarize(mass_tot75 = sum(mass))

# 75-100 cm depth sum
mass75_100_osu <- mass_osu %>%
  filter(hrzdep_t >=75 & hrzdep_t < 100) %>%
  group_by(dsp_pedon_id) %>%
  summarize(mass_tot100 = sum(mass))

# Join together into one dataframe and calculate min soil mass in each depth
mass_osu_all <- mass0_10_osu %>%
  left_join(mass10_30_osu, by="dsp_pedon_id") %>%
  left_join(mass30_50_osu, by="dsp_pedon_id") %>%
  left_join(mass50_75_osu, by="dsp_pedon_id") %>%
  left_join(mass75_100_osu, by="dsp_pedon_id")

mass_min_osu <- mass_osu_all %>%
  summarize(across(mass_tot10:mass_tot100, ~min(.x, na.rm=TRUE)))

# use minimum soil mass in each layer as reference mass, try manual calculation
# note - I don't know that using the minimum soil mass in each layer is the best way to choose a reference mass. I just needed something to work with.

# ESM calcs for OSU project with manual reference masses ----
fd_osu_manual <- read.csv(here("data_processed",  "Output_FD_osu_manual.csv"), sep = ";")
esm1_osu_manual <- read.csv(here("data_processed", "Output_ESM_osu_manual.csv"), sep = ";") %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE)
esm2_osu_manual <- read.csv(here("data_processed", "Output_ESM2_osu_manual.csv"), sep = ";") %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE)

min(esm1_osu_manual$Lower_cm)
# wahoo, deepest depth is now 90 cm :)

# calculate fixed depth SOC stocks for OSU project at intervals that match ESM calcs ----
# Calculate SOC in each horizon
soc_osu <- horizons(dice_osu) %>%
  mutate(hrzdepth = hrzdep_b - hrzdep_t,
         cf_mult = 1 - (coarse_frag_fill/100)) %>%
  mutate(soc_stock_hrz = soc_fill * bd_fill * hrzdepth * cf_mult)

# 0-10 cm depth sum
soc0_10_osu <- soc_osu %>%
  filter(hrzdep_t < 10) %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_tot10 = sum(soc_stock_hrz))

# 10-30 cm depth sum
soc10_30_osu <- soc_osu %>%
  filter(hrzdep_t >=10 & hrzdep_t < 30) %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_tot30 = sum(soc_stock_hrz))

# 30-50 cm depth sum
soc30_50_osu <- soc_osu %>%
  filter(hrzdep_t >=30 & hrzdep_t < 50) %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_tot50 = sum(soc_stock_hrz))

# 50-75 cm depth sum
soc50_75_osu <- soc_osu %>%
  filter(hrzdep_t >=50 & hrzdep_t < 75) %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_tot75 = sum(soc_stock_hrz))

# 75-100 cm depth sum
soc75_100_osu <- soc_osu %>%
  filter(hrzdep_t >=75 & hrzdep_t < 100) %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_tot100 = sum(soc_stock_hrz))

# join into one dataframe
soc_tot_osu <- soc0_10_osu %>%
  left_join(soc10_30_osu, by="dsp_pedon_id") %>%
  left_join(soc30_50_osu, by="dsp_pedon_id") %>%
  left_join(soc50_75_osu, by="dsp_pedon_id") %>%
  left_join(soc75_100_osu, by="dsp_pedon_id") %>%
  pivot_longer(cols = soc_tot10:soc_tot100,
               names_to = c(".value", "depth"),
               names_sep="_"
               ) %>%
  group_by(dsp_pedon_id) %>%
  mutate(depth_fd = ifelse(depth=="tot10", -10,
                        ifelse(depth=="tot30", -30,
                               ifelse(depth=="tot50", -50,
                                      ifelse(depth=="tot75", -75, -100))))) %>%
  mutate(layer=ifelse(depth=="tot10", 1,
                      ifelse(depth=="tot30", 2,
                             ifelse(depth=="tot50", 3,
                                    ifelse(depth=="tot75", 4, 5))))) %>%
  select(-depth) %>%
  unite("sample_id", c("dsp_pedon_id", "layer"), sep="-", remove=FALSE) %>%
  rename(soc_fd = soc) %>%
  ungroup()

# Join and compare OSU FD calcs to ESM calcs with manual reference masses ----

compare_osu_manual <- esm1_osu_manual %>%
  select(Campaign, Treatment, Point, sample_id,Layer, Lower_cm, SOC_stock_ESM) %>%
  rename(project = Campaign,
         label=Treatment,
         dsp_pedon_id = Point,
         depth_esm1 = Lower_cm,
         soc_esm1 = SOC_stock_ESM) %>%
  left_join(select(esm2_osu_manual, 
                   sample_id, Lower_cm, SOC_stock_ESM2), by="sample_id") %>%
  rename(soc_esm2 = SOC_stock_ESM2,
         depth_esm2 = Lower_cm) %>%
  left_join(select(soc_tot_osu, sample_id, soc_fd, depth_fd), by="sample_id") %>%
  pivot_longer(cols = depth_esm1:depth_fd,
               names_to = c(".value", "method"),
               names_sep="_"
               )

# Plot comparison of SOC stocks in each layer
# Labels for layers
layer_labels <- c("1" = "0-10 cm",
                  "2" = "10-30 cm",
                  "3" = "30-50 cm",
                  "4" = "50-75 cm",
                  "5" = "75-100 cm")

ggplot(compare_osu_manual %>% filter(dsp_pedon_id !="JoV2-3"), # filter out one sample that only has shallow values
       aes(x=label, y=soc, fill=method)) +
  geom_boxplot() +
  facet_wrap(~Layer, scales="free", labeller = labeller(Layer = layer_labels)) +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()

# Plot comparison of depths used
ggplot(compare_osu_manual%>% filter(dsp_pedon_id !="JoV2-3"), # filter out one sample that only has shallow values
       aes(x=label, y=depth, fill=method)) +
  geom_boxplot() +
  facet_wrap(~Layer, scales="free", labeller = labeller(Layer = layer_labels)) +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()

# Try to estimate horizon boundaries for a project ----
# Following along with this AQP tutorial: https://ncss-tech.github.io/AQP/aqp/estimate-ML-horizonation.html

# Subset for Kansas State project only
kansas <- soc_horizon %>%
  filter(project=="KansasState")

# Add generalized horizon labels (determined visually in "dsp4sh_prelim/code/02_profile_check.R")
# Sequence: A, Bt, Btk, Bk 
n_ks <- c('A', 'Bt', 'Btk', 'Bk') # generalized horizon label sequence
p_ks <- c('^A',
          '^Bt$|^Bt1|^Bt2|^Bt3|^Bt4',
          '^Btk',
          '^Bk')

# Run function to assign depths to generalized horizons
slice_vect <- seq(from = 0, to = 99, by = 1) # tell function what depths you want
dd <- NULL # need to run this or function will not work
depths_kansas_df <- genhz_depths(kansas, n_ks, p_ks)

# Use estimated horizon boundaries to calculate soil mass in each horizon ----
# Make vector of desired horizon depths
depths_kansas <- c(10,23,56,85,100)

# Calculated aggregated soil masses and save output
mass_agg_ks <- soil_mass_aggregate(kansas, depths_kansas)

# Compare Kansas ESM minimum mass calcs with FD calcs ----
# Read in calculations
fd_kansas_min <- read.csv(here("data_processed", "Output_FD_kansas_min.csv"), sep = ";")
esm1_kansas_min <- read.csv(here("data_processed", "Output_ESM_kansas_min.csv"), sep = ";") %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE)
esm2_kansas_min <- read.csv(here("data_processed", "Output_ESM2_kansas_min.csv"), sep = ";") %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE)

# Calculate fixed depth SOC stocks with depth increments that mirror ESM increments 
soc_agg_ks <- soc_stock_fd(kansas, depths_kansas)

# Join into one dataframe
compare_kansas_min <- esm1_kansas_min %>%
  select(Campaign, Treatment, Point, sample_id, Layer, Lower_cm, SOC_stock_ESM) %>%
  rename(project = Campaign,
         label=Treatment,
         dsp_pedon_id = Point,
         depth_esm1 = Lower_cm,
         soc_esm1 = SOC_stock_ESM) %>%
  left_join(select(esm2_kansas_min, 
                   sample_id, Lower_cm, SOC_stock_ESM2), by="sample_id") %>%
  rename(soc_esm2 = SOC_stock_ESM2,
         depth_esm2 = Lower_cm) %>%
  left_join(soc_agg_ks, by=c("sample_id", "dsp_pedon_id")) %>%
  select(-layer) %>%
  pivot_longer(cols = depth_esm1:depth_fd,
               names_to = c(".value", "method"),
               names_sep="_"
  ) %>%
  rename(soc_stock_calc = soc) %>%
  mutate(reference = "min")

# Plot comparison of SOC stocks in each layer
# Labels for layers
layer_labels_ks <- c("1" = "0-10 cm",
                  "2" = "10-23 cm",
                  "3" = "23-56 cm",
                  "4" = "56-85 cm",
                  "5" = "85-100 cm")

# Plot comparison of calculated SOC stocks
ggplot(compare_kansas_min,
       aes(x=label, y=soc_stock_calc, fill=method)) +
  geom_boxplot() +
  facet_wrap(~Layer, scales="free", labeller = labeller(Layer = layer_labels_ks)) +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()

# Plot comparison of depths used
ggplot(compare_kansas_min,
       aes(x=label, y=depth, fill=method)) +
  geom_boxplot() +
  facet_wrap(~Layer, scales="free", labeller = labeller(Layer = layer_labels_ks)) +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()

# Compare Kansas ESM maximum reference mass calcs with FD calcs ----
# This is sort of analogous to using the BAU as the reference case
# In almost every case, I imagine that Ref is going to have lower BD and therefore less mass in a fixed depth than BAU and SHM. Using BAU as the reference mass would require getting deeper in the Ref samples...
# can test this by using the max values as ESM input
# Read in calculations
esm1_kansas_max <- read.csv(here("data_processed", "Output_ESM_kansas_max.csv"), sep = ";") %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE)
esm2_kansas_max <- read.csv(here("data_processed", "Output_ESM2_kansas_max.csv"), sep = ";") %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE)

# Join into one dataframe
compare_kansas_max <- esm1_kansas_max %>%
  select(Campaign, Treatment, Point, sample_id, Layer, Lower_cm, SOC_stock_ESM) %>%
  rename(project = Campaign,
         label=Treatment,
         dsp_pedon_id = Point,
         depth_esm1 = Lower_cm,
         soc_esm1 = SOC_stock_ESM) %>%
  left_join(select(esm2_kansas_max, 
                   sample_id, Lower_cm, SOC_stock_ESM2), by="sample_id") %>%
  rename(soc_esm2 = SOC_stock_ESM2,
         depth_esm2 = Lower_cm) %>%
  left_join(soc_agg_ks, by=c("sample_id", "dsp_pedon_id")) %>%
  select(-layer) %>%
  pivot_longer(cols = depth_esm1:depth_fd,
               names_to = c(".value", "method"),
               names_sep="_"
  ) %>%
  rename(soc_stock_calc = soc) %>%
  mutate(reference = "max")

# Plot comparison of calculated SOC stocks
ggplot(compare_kansas_max,
       aes(x=label, y=soc_stock_calc, fill=method)) +
  geom_boxplot() +
  facet_wrap(~Layer, scales="free", labeller = labeller(Layer = layer_labels_ks)) +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()

# Plot comparison of depths used
ggplot(compare_kansas_max,
       aes(x=label, y=depth, fill=method)) +
  geom_boxplot() +
  facet_wrap(~Layer, scales="free", labeller = labeller(Layer = layer_labels_ks)) +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()

# How was function even able to calculate ESM masses for soil depths that are below what was sampled? I get that the cubic spline can probably do that for the ESM2 calcs, but what about classical ESM?

# Compare Kansas ESM average calcs with FD calcs ----
# Read in calculations
esm1_kansas_avg <- read.csv(here("data_processed", "Output_ESM_kansas_avg.csv"), sep = ";") %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE)
esm2_kansas_avg <- read.csv(here("data_processed", "Output_ESM2_kansas_avg.csv"), sep = ";") %>%
  unite("sample_id", Point:Layer, sep="-", remove=FALSE)

# Join into one dataframe
compare_kansas_avg <- esm1_kansas_avg %>%
  select(Campaign, Treatment, Point, sample_id, Layer, Lower_cm, SOC_stock_ESM) %>%
  rename(project = Campaign,
         label=Treatment,
         dsp_pedon_id = Point,
         depth_esm1 = Lower_cm,
         soc_esm1 = SOC_stock_ESM) %>%
  left_join(select(esm2_kansas_avg, 
                   sample_id, Lower_cm, SOC_stock_ESM2), by="sample_id") %>%
  rename(soc_esm2 = SOC_stock_ESM2,
         depth_esm2 = Lower_cm) %>%
  left_join(soc_agg_ks, by=c("sample_id", "dsp_pedon_id")) %>%
  select(-layer) %>%
  pivot_longer(cols = depth_esm1:depth_fd,
               names_to = c(".value", "method"),
               names_sep="_"
  ) %>%
  rename(soc_stock_calc = soc) %>%
  mutate(reference = "avg")

# Plot comparison of calculated SOC stocks
ggplot(compare_kansas_avg,
       aes(x=label, y=soc_stock_calc, fill=method)) +
  geom_boxplot() +
  facet_wrap(~Layer, scales="free", labeller = labeller(Layer = layer_labels_ks)) +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()

# Plot comparison of depths used
ggplot(compare_kansas_avg,
       aes(x=label, y=depth, fill=method)) +
  geom_boxplot() +
  facet_wrap(~Layer, scales="free", labeller = labeller(Layer = layer_labels_ks)) +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()

# Put all Kansas ESM methods together (min, max, average) ----
kansas_esm_all <- rbind(compare_kansas_min, compare_kansas_max, compare_kansas_avg)

# Compare 0-10 cm SOC stocks via all reference methods
ggplot(kansas_esm_all %>% filter(Layer==1),
       aes(x=label, y=soc_stock_calc, fill=method)) +
  geom_boxplot() +
  facet_wrap(~reference) +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()