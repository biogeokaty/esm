# ESM Methods Comparison Analysis

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
library(lmerTest)
library(multcomp)
library(ggeffects)
library(tidymodels)
library(multilevelmod)
library(broom.mixed)
library(flextable)
library(tidyverse)

# Data
esm_all <- read.csv(here("data_processed", "esm_all.csv"))
soc_horizon <- read.csv(here("data_processed","04_soc_stock_horizon.csv"))

# Useful subsets of data
esm_standard <- esm_all %>%
  filter(depth_increments == "standard", ref_data=="indv_project")

esm_gen <- esm_all %>%
  filter(depth_increments == "genetic_hrz", ref_data=="indv_project")

esm_same <- esm_all %>%
  filter(ref_data=="all_data")

# Make list of summary functions - useful for generating summary tables
sd_cv <- list(
  sd = ~round(sd(.x, na.rm = TRUE), 2),
  cv = ~round(((sd(.x, na.rm=TRUE) / mean(.x, na.rm=TRUE))* 100), 2)
)

# list of projects - will need for purrr later
projects <- esm_all %>% 
  distinct(project) %>% 
  pull()

# depth layer labels
std_depth_labels <- c(
  "1" = "0-10 cm",
  "2" = "10-30 cm",
  "3" = "30-50 cm",
  "4" = "50-75 cm",
  "5" = "75-100 cm"
)

# Exploratory figures  - how does ESM choice influence total calculated SOC stocks? ----
# Compare calculated 0-10 cm SOC stocks in all projects, gridded by method and colored by reference mass summary statistic
# This will be the same for genetic horizon and standard depth increment data, so we can choose either to plot 

# Try looking at all projects together
ggplot(esm_standard %>% filter(layer==1),
       aes(x=label, y=soc, fill=ref_stat)) +
  geom_boxplot() +
  facet_wrap(~method) +
  labs(y="SOC Stock (Mg/ha) in 0-10 cm") +
  scale_fill_paletteer_d("nationalparkcolors::Arches", name="Reference Mass") +
  theme_classic()
# Hard to read/interpret this since the projects have a lot of variability in SOC - try to plot them individually

# Plot 0-10cm stocks for all projects
proj_surf_method_plots <- purrr::map(.x = projects,
                                     .f = ~{
                                       esm_standard %>% 
                                         filter(layer == 1, project==.x) %>%
                                         ggplot(aes(x=label, y=soc, fill=ref_stat)) +
                                         geom_boxplot() +
                                         facet_wrap(~method) +
                                         labs(x="Management",
                                              y="SOC Stock (Mg/ha) in 0-10 cm",
                                              title=glue::glue(.x, " - 0-10 cm SOC Stocks")) +
                                         scale_fill_paletteer_d("nationalparkcolors::Arches", name="Reference Mass") +
                                         theme_classic() +
                                         theme(plot.title=element_text(hjust=0.5))
                                       
                                       ggsave(here("figs", glue::glue("surface_esm", .x, ".png")), 
                                              width=10, height=7, units="in", dpi=400)
                                     })

# Calculate 0-30 cm and 0-100 cm totals for standard depth increment data
esm_standard_totals <- esm_standard %>%
  group_by(project, depth_increments, method_long, method, ref_stat, label, dsp_pedon_id) %>%
  summarize(soc_0to10 = sum(soc[apparent_depth=="0-10 cm"]),
            soc_0to30 = sum(soc[apparent_depth=="0-10 cm" | apparent_depth=="10-30 cm"]),
            soc_0to100 = sum(soc))

esm_standard_totals_longer <- esm_standard_totals %>%
  pivot_longer(cols=soc_0to10:soc_0to100,
               names_to = c(".value", "depth"),
               names_sep="_") %>%
  mutate(label=factor(label, levels=c("BAU", "SHM", "Ref")),
         depth=factor(depth, levels=c("0to10", "0to30", "0to100")))

# Plot total SOC stocks in different depth increments for different methods
# Example - Illinois
esm_totals_ill <- esm_standard_totals_longer %>%
  filter(project=="Illinois")

ggplot(esm_totals_ill, aes(x=label, y=soc, fill=ref_stat)) +
  geom_boxplot() +
  facet_grid(depth ~ method, scales="free_y") +
  scale_fill_paletteer_d("nationalparkcolors::Arches", name="Reference Mass") +
  labs(title="Comparison of ESM Methods for \nCalcuating SOC Stock Totals - Illinois",
       y = "SOC Stock (Mg/ha)",
       x ="Management") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5))

# Make with purrr for all projects
stock_totals_plots <- purrr::map(.x = projects,
                                 .f = ~{
                                   esm_standard_totals_longer %>%
                                     filter(project==.x) %>%
                                     ggplot(aes(x=label, y=soc, fill=ref_stat)) +
                                     geom_boxplot() +
                                     facet_grid(depth ~ method, scales="free_y") +
                                     scale_fill_paletteer_d("nationalparkcolors::Arches", name="Reference Mass") +
                                     labs(title=glue::glue("Comparison of ESM Methods for \nCalcuating SOC Stock Totals - ", .x),
                                          y = "SOC Stock (Mg/ha)",
                                          x ="Management") +
                                     theme_classic() +
                                     theme(plot.title=element_text(hjust=0.5))
                                   
                                   ggsave(here("figs", glue::glue("stock_totals_", .x, ".png")), 
                                          width=10, height=7, units="in", dpi=400)
                                 })

# What is the appropriate statistic to test?
# Mixed linear model seems appropriate....
# Fixed effect is calculation method, random effect is project and label (to account for variation in SOC stocks by project and management)

# For 0-10 cm stocks
esm_std_0to10_lmer <- lmer(soc_0to10 ~ method_long + (1|project/label), data = esm_standard_totals)
summary(esm_std_0to10_lmer)
plot(esm_std_0to10_lmer)
qqnorm(resid(esm_std_0to10_lmer))
qqline(resid(esm_std_0to10_lmer))
# QQPlot indicates we are missing some source of variation...

# Is method significant overall?
drop1(lmer(soc_0to10 ~ method_long + (1|project/label), data = esm_standard_totals, REML=FALSE), test = "Chisq")
# yes

# For 0-30 cm stocks
esm_std_0to30_lmer <- lmer(soc_0to30 ~ method_long + (1|project/label), data = esm_standard_totals)
summary(esm_std_0to30_lmer)
plot(esm_std_0to30_lmer)
qqnorm(resid(esm_std_0to30_lmer))
qqline(resid(esm_std_0to30_lmer))

# Is method significant overall?
drop1(lmer(soc_0to30 ~ method_long + (1|project/label), data = esm_standard_totals, REML=FALSE), test = "Chisq")
# yes

esm_std_0to100_lmer <- lmer(soc_0to100 ~ method_long + (1|project/label), data = esm_standard_totals)
summary(esm_std_0to100_lmer)
plot(esm_std_0to100_lmer)
qqnorm(resid(esm_std_0to100_lmer))
qqline(resid(esm_std_0to100_lmer))

# Is method significant overall?
drop1(lmer(soc_0to100 ~ method_long + (1|project/label), data = esm_standard_totals, REML=FALSE), test = "Chisq")
# yes

# Mixed linear models for ESM data only - test method (1 vs 2) and reference stat (min, mean, max)
esm_standard_totals_no_fd <- esm_standard_totals %>%
  filter(method!="fd")

# Fit model
esm_std_0to10_no_fd_lmer <- lmer(soc_0to10 ~ method + ref_stat + (1|project/label), data = esm_standard_totals_no_fd)
summary(esm_std_0to10_no_fd_lmer)

# test significance of fixed effects
drop1(lmer(soc_0to10 ~ method + ref_stat + (1|project/label), data = esm_standard_totals_no_fd, REML=FALSE), test = "Chisq")
# no difference between ESM1 and ESM2, but significant effect of reference stat

# Repeat for other depths
# 0-30 cm depth
esm_std_0to30_no_fd_lmer <- lmer(soc_0to30 ~ method + ref_stat + (1|project/label), data = esm_standard_totals_no_fd)
summary(esm_std_0to30_no_fd_lmer)

# test significance of fixed effects
drop1(lmer(soc_0to30 ~ method + ref_stat + (1|project/label), data = esm_standard_totals_no_fd, REML=FALSE), test = "Chisq")
# method is nearly significant, significant effect of reference stat

# 0-100 cm depth
esm_std_0to100_no_fd_lmer <- lmer(soc_0to100 ~ method + ref_stat + (1|project/label), data = esm_standard_totals_no_fd)
summary(esm_std_0to100_no_fd_lmer)

# test significance of fixed effects
drop1(lmer(soc_0to100 ~ method + ref_stat + (1|project/label), data = esm_standard_totals_no_fd, REML=FALSE), test = "Chisq")
# significant effect of reference stat only

# How does ESM method choice influence total stock error? ----
esm2_std_totals <- esm_standard_totals_longer %>%
  filter(method=="esm2")

esm2_std_totals_error <- esm2_std_totals %>%
  group_by(depth, ref_stat) %>%
  summarize(across(soc, sd_cv)) 
flextable(esm2_std_totals_error)

# How does ESM method choice influence depths used to calculate SOC stocks ----
  
# Plot comparison of depths used - this gives a good sense of how reasonable the data is IMO
# try for one project
ggplot(esm_standard %>% filter(project=="Illinois", method=="esm2"), aes(x=ref_stat, y=depth, fill=ref_stat)) +
  geom_boxplot() +
  facet_wrap(~layer, scales="free_y", labeller=labeller(layer=std_depth_labels)) +
  labs(x="Reference Mass",
       y="Actual bottom depth used for ESM calculations", 
       title="Influence of Reference Mass on Soil Depths \nUsed for ESM Cubic Spline Calculations - Illinois") +
  scale_fill_paletteer_d("nationalparkcolors::Arches") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5),
        legend.position="none")

# also look at what happens if we use the same reference mass for all projects
ggplot(esm_all %>% filter(depth_increments=="standard", ref_data=="all_data", 
                          project=="Illinois", method=="esm2"), aes(x=ref_stat, y=depth, fill=ref_stat)) +
  geom_boxplot() +
  facet_wrap(~layer, scales="free_y", labeller=labeller(layer=std_depth_labels)) +
  labs(x="Reference Mass",
       y="Actual bottom depth used for ESM calculations", 
       title="Influence of Reference Mass on Soil Depths \nUsed for ESM Cubic Spline Calculations - Illinois") +
  scale_fill_paletteer_d("nationalparkcolors::Arches") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5),
        legend.position="none")

# can we plot these together?
ggplot(esm_all %>% filter(depth_increments=="standard", project=="Illinois", method=="esm2"), 
       aes(x=ref_stat, y=depth, fill=ref_data)) +
  geom_boxplot() +
  facet_wrap(~layer, scales="free_y", labeller=labeller(layer=std_depth_labels)) +
  labs(x="Reference Summary Statistic",
       y="Actual bottom depth used for ESM calculations", 
       title="Influence of Reference Mass on Soil Depths \nUsed for ESM Cubic Spline Calculations - Illinois") +
  scale_fill_paletteer_d("nationalparkcolors::Arches", 
                         name="Reference Data Source", labels=c("Entire DSP4SH dataset", "Individual project")) +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5))

ggplot(esm_all %>% filter(depth_increments=="standard", project=="OregonState", method=="esm2"), 
       aes(x=ref_stat, y=depth, fill=ref_data)) +
  geom_boxplot() +
  facet_wrap(~layer, scales="free_y", labeller=labeller(layer=std_depth_labels)) +
  labs(x="Reference Summary Statistic",
       y="Actual bottom depth used for ESM calculations", 
       title="Influence of Reference Mass on Soil Depths \nUsed for ESM Cubic Spline Calculations - Oregon State") +
  scale_fill_paletteer_d("nationalparkcolors::Arches", 
                         name="Reference Data Source", labels=c("Entire DSP4SH dataset", "Individual project")) +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5))

# Plot for all projects
purrr::map(.x = projects,
           .f = ~{
             cubic <- esm_all %>% 
               filter(depth_increments=="standard", method=="esm2", project==.x) %>%
               ggplot(aes(x=ref_stat, y=depth, fill=ref_data)) +
               geom_boxplot() +
               facet_wrap(~layer, scales="free_y", labeller=labeller(layer=std_depth_labels)) +
               labs(x="Reference Summary Statistic",
                    y="Actual bottom depth used for ESM calculations", 
                    title="Cubic spline ESM") +
               scale_fill_paletteer_d("nationalparkcolors::Arches", 
                                      name="Reference Data Source", labels=c("Entire DSP4SH dataset", "Individual project")) +
               theme_classic() +
               theme(plot.title=element_text(hjust=0.5))
             
             classic <- esm_all %>% 
               filter(depth_increments=="standard", method=="esm1", project==.x) %>%
               ggplot(aes(x=ref_stat, y=depth, fill=ref_data)) +
               geom_boxplot() +
               facet_wrap(~layer, scales="free_y", labeller=labeller(layer=std_depth_labels)) +
               labs(x="Reference Summary Statistic",
                    y="Actual bottom depth used for ESM calculations", 
                    title="Classic (non-modeled) ESM") +
               scale_fill_paletteer_d("nationalparkcolors::Arches", 
                                      name="Reference Data Source", labels=c("Entire DSP4SH dataset", "Individual project")) +
               theme_classic() +
               theme(plot.title=element_text(hjust=0.5))
             
             plot_row <- plot_grid(classic + theme(legend.position="none"), 
                                   cubic, 
                                   nrow=1, labels=c("A", "B"), rel_widths=c(1, 1.3))
             
             title <- ggdraw() +
                     draw_label(glue::glue("Influence of Reference Mass Choice on Soil Depths Used for ESM Calculations - ", .x))
             
             plot_grid(title, plot_row, ncol=1, rel_heights = c(0.1, 1))
             
             ggsave(here("figs", glue::glue("esm_depths_", .x, ".png")), 
                    width=14, height=10, units="in", dpi=400)

           })

# there are some very funky depths!
# pull out which pedons, and if the depths are off when using reference masses for individual projects or the entire dataset
deep <- esm_all %>%
  filter(depth < -220) %>%
  distinct(dsp_pedon_id, ref_data, ref_stat) 

# funky depths only happen when using the same reference mass for the entire dataset, and almost always when using the maximum soil mass as a reference

# for one project - how do standard depths vs genetic horizons compare?
# make plots
std_depth_plots <- purrr::map(.x = projects,
                               .f = ~{
                                 esm_standard %>% 
                                   filter(method == "esm2", project==.x) %>%
                                   ggplot(aes(x=label, y=soc, fill=ref_stat)) +
                                   geom_boxplot() +
                                   facet_wrap(~apparent_depth, scales="free") +
                                   labs(x="Management",
                                        y="SOC Stock (Mg/ha)",
                                        title=glue::glue(.x, " - Standard Depth Increments")) +
                                   scale_fill_paletteer_d("nationalparkcolors::Arches", name="Reference Mass") +
                                   theme_classic() +
                                   theme(plot.title=element_text(hjust=0.5))
                                 
                                 ggsave(here("figs", glue::glue("esm_standard_", .x, ".png")), 
                                        width=10, height=7, units="in", dpi=400)
                               })

# Compare within-layer error for standard vs genetic horizon depth increments
# calculate CV for each layer within each project/management group?
depth_comp_error <- esm_all %>%
  filter(ref_data=="indv_project") %>%
  group_by(project, label, method_long, depth_increments, layer) %>%
  summarize(layer_cv = round(((sd(soc, na.rm=TRUE) / mean(soc, na.rm=TRUE))* 100), 2))

# Plot CVs
ggplot(depth_comp_error, aes(y=layer_cv, x=method_long, fill=depth_increments)) +
  geom_boxplot() +
  facet_grid(project~layer, scales="free_y") +
  labs(x="SOC stock calculation method", y="Layer Coefficient of Variation") +
  scale_fill_paletteer_d("nationalparkcolors::Arches", name="Depth Increments") +
  theme_classic()
# This figure is small and hard to interpret but overall I am pretty sure it shows that the within-layer CVs are pretty much the same whether SOC stocks are calculated on a standard depth increments basis or a genetic horizon basis

# try to pull out just one project...
ggplot(depth_comp_error %>% filter(project=="Illinois"), aes(y=layer_cv, x=method_long, fill=depth_increments)) +
  geom_boxplot() +
  facet_wrap(~layer, scales="free_y") +
  labs(x="SOC stock calculation method", y="Layer Coefficient of Variation", title="Illinois") +
  scale_fill_paletteer_d("nationalparkcolors::Arches", name="Depth Increments") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
# perhaps there is slightly greater error in the ESM2 method vs ESM1 method, especially in greater depths - because of more extrapolation?

# can figure out if SOC stock calculation method contributes more error than depth increments with an lme
depth_comp_lmer <- lmer(layer_cv ~ method_long + depth_increments + (1|project/label/layer), data = depth_comp_error)
summary(depth_comp_lmer)
plot(depth_comp_lmer)
qqnorm(resid(depth_comp_lmer))
qqline(resid(depth_comp_lmer))

# Is method or depth increment significant overall?
drop1(lmer(layer_cv ~ method_long + depth_increments + (1|project/label/layer), data = depth_comp_error, REML=FALSE), test = "Chisq")
# Method significantly influences within-layer error, but depth increment choice does not

# should be able to figure out which method results in least within-layer error????
summary(glht(depth_comp_lmer, linfct = mcp(method_long = 'Tukey')))

# Plot predictions of within-layer error for each method
depth_comp_pred <- ggpredict(depth_comp_lmer, terms = c("method_long"))

ggplot(depth_comp_pred , aes(x=x, y=predicted)) +
  geom_point() +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high)) +
  labs(x="SOC stock calculation method", y="Predicted within-layer CV") +
  theme_classic()

# Does choice in depth increments alter total stocks calculated? ----
# Calculate 0 to 10 cm and 0 to 100 cm SOC stocks for genetic horizon depth increment data
# the 0-10 cm SOC stocks absolutely should be the same between methods so this is really just a check

esm_gen_totals <- esm_gen %>%
  group_by(project, depth_increments, method_long, method, ref_stat, label, dsp_pedon_id) %>%
  summarize(soc_0to10 = sum(soc[apparent_depth=="0-10 cm"]),
            soc_0to100 = sum(soc))

esm_gen_totals_longer <- esm_gen_totals %>%
  pivot_longer(cols=soc_0to10:soc_0to100,
               names_to = c(".value", "depth"),
               names_sep="_") %>%
  mutate(label=factor(label, levels=c("BAU", "SHM", "Ref")),
         depth=factor(depth, levels=c("0to10", "0to100")))

esm_gen_std_totals <- rbind(esm_standard_totals_longer, esm_gen_totals_longer)

# plot for a project
ggplot(esm_gen_std_totals %>% filter(project == "Illinois", depth!="0to30"), aes(x=method_long, y=soc, fill=depth_increments)) +
  geom_boxplot() +
  facet_wrap(~depth, scales="free_y") +
  labs(x="SOC stock calculation method",
       y="SOC Stock (Mg/ha)",
       title="Influence of depth increment selection on \nSOC stock calculation - Illinois") +
  scale_fill_paletteer_d("nationalparkcolors::Arches", name="Depth Increments") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5),
        axis.text.x=element_text(angle=45, hjust=1))

# Depth increment selection does NOT alter calculated SOC stocks - this is expected (it's the same data, after all!), but is just a good check that we aren't creating crazy errors by re-jiggering the data into standard depth increments. Huzzah!  

# try for other projects
purrr::map(.x = projects, .f= ~{
  ggplot(esm_gen_std_totals %>% filter(project == .x, depth!="0to30"), 
         aes(x=method_long, y=soc, fill=depth_increments)) +
    geom_boxplot() +
    facet_wrap(~depth, scales="free_y") +
    labs(x="SOC stock calculation method",
         y="SOC Stock (Mg/ha)",
         title=glue::glue("Influence of depth increment selection on \nSOC stock calculation - ", .x)) +
    scale_fill_paletteer_d("nationalparkcolors::Arches", name="Depth Increments") +
    theme_classic() +
    theme(plot.title=element_text(hjust=0.5),
          axis.text.x=element_text(angle=45, hjust=1))
})

# How does ESM method choice influence between-group differences in SOC stocks? ----
# Focus on standard depth increments here because we know there is no significant difference in stocks based on depth increment method selection, and standard depth increments are easier to compare between projects

# Total stocks: 0-10 cm
# for each different method, do we see significant effect of project?
# fixed depth
fd <- esm_standard_totals %>%
  filter(method=="fd")

fd_0to10_lmer <- lmer(soc_0to10 ~ label + (1|project), data = esm_standard_totals)
summary(fd_0to10_lmer)
drop1(lmer(soc_0to10 ~ label + (1|project), data = esm_standard_totals), REML=FALSE, test="Chisq")
# overall model is significant - management does significantly influence 0-10cm SOC stocks (FD) once between-project variation is accounted for

# which groups are significantly different
summary(glht(fd_0to10_lmer, linfct = mcp(label = 'Tukey')))
# All groups are significantly different

# Plot predictions of SOC stocks for each method
fd_0to10_pred <- ggpredict(fd_0to10_lmer, terms = c("label"))

ggplot(fd_0to10_pred, aes(x=x, y=predicted)) +
  geom_point() +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high)) +
  labs(x="Management", y="Predicted 0-10 cm SOC stock (Mg/ha)") +
  theme_classic()

# Run linear mixed model and extract significance
esm_group_comp <- esm_standard_totals_longer %>%
  group_by(method_long, depth) %>%
  nest() %>%
  mutate(lmer = map(data, ~lmer(soc ~ label + (1|project), data = .x)),
  drop1 = map(data, .f = ~{
    drop1(lmer(soc ~ label + (1|project), REML=FALSE, data = .x))
    }),
  tidy_lmer = map(lmer, broom.mixed::tidy),
  tidy_drop = map(drop1, broom.mixed::tidy)) %>%
  select(method_long, depth, tidy_lmer, tidy_drop) %>%
  unnest(cols=c(tidy_lmer, tidy_drop), names_sep="_") %>%
  mutate(sig = case_when(tidy_drop_p.value < 0.05 ~ "significant",
                         tidy_drop_p.value > 0.05 ~ "not_significant"))

total_not_sig <- esm_group_comp %>%
  filter(sig=="not_significant") %>%
  distinct(method_long, depth, sig)

# The only method that doesn't detect significant between-group differences is ESM2 (any reference soil choice) for the 0 to 100cm total stocks - that seems reasonable

# Layer-specific stocks
# Run linear mixed model and extract significance
esm_layer_group_comp <- esm_standard %>%
  group_by(method_long, layer) %>%
  nest() %>%
  mutate(lmer = map(data, ~lmer(soc ~ label + (1|project), data = .x)),
         drop1 = map(data, .f = ~{
           drop1(lmer(soc ~ label + (1|project), REML=FALSE, data = .x))
         }),
         tidy_lmer = map(lmer, broom.mixed::tidy),
         tidy_drop = map(drop1, broom.mixed::tidy)) %>%
  select(method_long, layer, tidy_lmer, tidy_drop) %>%
  unnest(cols=c(tidy_lmer, tidy_drop), names_sep="_") %>%
  mutate(sig = case_when(tidy_drop_p.value < 0.05 ~ "significant",
                         tidy_drop_p.value > 0.05 ~ "not_significant"))

layer_not_sig <- esm_layer_group_comp %>%
  filter(sig=="not_significant") %>%
  distinct(method_long, layer, sig) %>%
  arrange(layer)
# all methods find significant effects of management in 0-10cm and 10-30 cm layers
# ESM2 with max reference mass finds no sig in 30-50 cm
# most calculation approaches find no significant management effect in 50-75cm (exceptions here: fixed depth and ESM1 min)
# no calculation approaches find a significant management effect 75-100 cm layers