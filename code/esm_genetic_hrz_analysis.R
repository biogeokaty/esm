# ESM Analysis - Genetic Horizon Increments

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

# Nest the data ----
esm_gen_nested <- esm_gen %>%
  group_by(project, ref_stat, depth_increments) %>%
  nest()

# list of projects
projects <- esm_gen %>% 
  distinct(project) %>% 
  pull()

# Preliminary figures ----

# Compare calculated 0-10 cm SOC stocks in all projects, gridded by method and colored by reference mass summary statistic
ggplot(esm_gen %>% filter(layer==1),
       aes(x=label, y=soc, fill=ref_stat)) +
  geom_boxplot() +
  facet_wrap(~method) +
  labs(y="SOC Stock (Mg/ha) in 0-10 cm") +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()

proj_surf_method_plots <- purrr::map(.x = projects,
                                     .f = ~{
                                       esm_gen %>% 
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

# try this for all projects individually

# Just look at ESM2 (cubic spline) method, all depths
ggplot(esm_gen %>% filter(method=="esm2"),
       aes(x=label, y=soc, fill=ref_stat)) +
  geom_boxplot() +
  facet_wrap(~layer, scales="free") +
  labs(y="SOC Stock (Mg/ha)") +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()
# this is really hard to look at with all the projects together

# make individual project plots using map() function
proj_depth_plots <- purrr::map(.x = projects,
                 .f = ~{
                   esm_gen %>% 
                     filter(method == "esm2", project==.x) %>%
                     ggplot(aes(x=label, y=soc, fill=ref_stat)) +
                     geom_boxplot() +
                     facet_wrap(~apparent_depth, scales="free") +
                     labs(x="Management",
                          y="SOC Stock (Mg/ha)",
                          title=glue::glue(.x, " - Depth Increments by Genetic Horizon")) +
                     scale_fill_paletteer_d("nationalparkcolors::Arches", name="Reference Mass") +
                     theme_classic() +
                     theme(plot.title=element_text(hjust=0.5))
                   
                   ggsave(here("figs", glue::glue("esm_", .x, ".png")), 
                          width=10, height=7, units="in", dpi=400)
                 })

plot_grid(plotlist = proj_depth_plots)
ggsave(here("figs", "esm2_genetic_depth_plots.png"), width=16, height=14, units="in", dpi=400)

# Initial statistical questions ----
# Is there a difference between ESM1 and ESM2? ESM and FD?
summary(lmer(soc ~ method + (1|project) + (1|layer) + (1|label) + (1|ref_stat), data = esm_gen))
full <- lmer(soc ~ method + (1|project) + (1|layer) + (1|label) + (1|ref_stat), data = esm_gen, REML = FALSE)
reduced <- lmer(soc ~ (1|project) + (1|layer) + (1|label) + (1|ref_stat), data = esm_gen, REML = FALSE)
anova(full, reduced, test="Chisq")
# significant effect of method

full2 <- lmer(soc ~ method + ref_stat + (1|project) + (1|layer) + (1|label), data = esm_gen, REML = FALSE)
reduced2 <- lmer(soc ~ (1|project) + (1|layer) + (1|label), data = esm_gen, REML = FALSE)
anova(full2, reduced2, test="Chisq")
# significant effect of ref stat

# maybe should run on subsets?? like just the 0-10cm??