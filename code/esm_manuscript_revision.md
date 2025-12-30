esm_manuscript_revision
================
Katy Dynarski
2025-12-30

# Overview

## Calculation method comparison

I calculated SOC stocks in 0-5 cm, 5-10 cm, 10-30 cm, and 30-60 cm depth
increments via fixed depth and ESM methods. For the ESM calculations, I
tested four different reference mass options:

- Minimum soil mass in an individual DSP4SH project (ESM project, min)
- Mean soil mass in an individual DSP4SH project (ESM project, mean)
- Minimum soil mass in the reference treatment in an individual DSP4SH
  project (ESM treatment, min)
- Mean soil mass in the reference treatment in an individual DSP4SH
  project (ESM treatment, mean)

I wanted to test just using data from the reference treatment to
calculate the reference mass because a number of studies calculated
their reference mass as either the mean or minimum soil mass in the
treatment expected to have the lowest bulk density.

I also re-calculated bulk density after reading through the SimpleESM
documentation again. We did not originally account for coarse fragment
content, so I re-calculated bulk density to correct for coarse fragment
content. It looks like this resolved some of the weird numbers we were
getting for soils high in coarse frags (e.g. the Oregon State project).

# Climate and soils data for each project

``` r
# need to rename OSU projects in "project" dataframe and calculate mean MAT/MAP for all projects
project_renamed <- project %>%
  mutate(project = case_when(
    project=="OregonState" & soil=="Jory" ~ "OregonStateJory",
    project=="OregonState" & soil=="Woodburn" ~ "OregonStateWoodburn", 
    .default = project))

project_avg_climate <- project_renamed %>%
  group_by(project) %>%
  dplyr::summarize(avg_mat = round(mean(mat, na.rm=TRUE), 1),
         avg_ppt = round(mean(map, na.rm=TRUE), 0)) %>%
  arrange(factor(project, levels=project_plotting_order))

# write_csv(project_avg_climate, here("figs", "revision_figs", "tablesupp1_site_clim.csv"))
```

# Number of pedons and sites

``` r
# how many projects
n_distinct(esm$project)
```

    ## [1] 9

``` r
# how many pedons
n_distinct(esm$dsp_pedon_id)
```

    ## [1] 205

``` r
# how many sites
# have to use some regex to extract the site name from the pedon ID since I'm missing the needed column
sites <- esm %>%
  select(project, label, dsp_pedon_id) %>%
  distinct(project, label, dsp_pedon_id) %>%
  mutate(dsp_site_id = str_remove(dsp_pedon_id, "-\\d*\\w*$")) %>%
  group_by(project, label, dsp_site_id)

# how many pedons per site?
sites %>%
  dplyr::summarize(n_pedons_per_site = n())
```

    ## `summarise()` has grouped output by 'project', 'label'. You can override using
    ## the `.groups` argument.

    ## # A tibble: 69 × 4
    ## # Groups:   project, label [27]
    ##    project  label dsp_site_id     n_pedons_per_site
    ##    <fct>    <fct> <chr>                       <int>
    ##  1 UConn    BAU   PaT2                            1
    ##  2 UConn    SHM   PaH2                            1
    ##  3 UConn    SHM   PaL1                            1
    ##  4 UConn    SHM   PaM2                            1
    ##  5 UConn    Ref   PaU1                            1
    ##  6 UConn    Ref   PaU2                            1
    ##  7 Illinois BAU   IL-115-SW-ETAL                  3
    ##  8 Illinois BAU   IL-115-SW-WDFLD                 3
    ##  9 Illinois BAU   IL-147-SW-CRP                   3
    ## 10 Illinois BAU   IL-147-SW-TS                    3
    ## # ℹ 59 more rows

``` r
# how many total sites  
n_distinct(sites$dsp_site_id)
```

    ## [1] 69

``` r
# how many total pedons
n_distinct(sites$dsp_pedon_id)
```

    ## [1] 205

``` r
# how many total pedons per treatment
esm %>%
  select(project, label, dsp_pedon_id) %>%
  distinct(project, label, dsp_pedon_id) %>%
  group_by(label) %>%
  dplyr::summarize(n_pedons_per_treatment = n())
```

    ## # A tibble: 3 × 2
    ##   label n_pedons_per_treatment
    ##   <fct>                  <int>
    ## 1 BAU                       76
    ## 2 SHM                       77
    ## 3 Ref                       52

``` r
# table for supplement with number of sites and pedons per treatment in each project
site_table <- sites %>%
  distinct(project, label, dsp_site_id) %>%
  group_by(project, label) %>%
  dplyr::summarize(n_sites_per_treatment = n())
```

    ## `summarise()` has grouped output by 'project'. You can override using the
    ## `.groups` argument.

``` r
pedon_table <- sites %>%
  ungroup() %>%
  group_by(project, label) %>%
  dplyr::summarize(n_pedons_per_trt = n())
```

    ## `summarise()` has grouped output by 'project'. You can override using the
    ## `.groups` argument.

``` r
site_pedon_table <- site_table %>%
  left_join(pedon_table, by=c("project", "label")) %>%
  arrange(factor(project, levels=project_plotting_order), label)

# write_csv(site_pedon_table, here("figs", "revision_figs", "tablesupp2_site_pedon_numbers.csv"))
```

# Range in bulk density and reference masses

Management can influence soil bulk density, which results in different
masses of soil within the same depth increment between treatments. We
can check this by plotting soil bulk density under the different
management treatments in each project.

Plot bulk density:

``` r
# Plot bulk density for each project/treatment
ggplot(bd_depth, aes(x=depth, y=depth_wt_bd, fill=label)) +
  geom_boxplot(fatten=1.5, lwd=0.3) +
  facet_wrap(~project, scales="free_y", labeller=labeller(project=project_labels_esm)) +
  scale_fill_manual(values=c("#FED789FF","#72874EFF","#476F84FF"),
                     breaks=c("BAU", "SHM", "Ref"), 
                    name="Management") +
  labs(x="Depth", y=expression("Bulk density"~(g ~ cm^-3))) +
  theme_classic() +
  theme(axis.text.x=element_text(hjust=1, angle=45))
```

    ## Warning: Removed 4 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

![](esm_manuscript_revision_files/figure-gfm/fig%20s2%20plot%20bulk%20density%20by%20depth%20increment-1.png)<!-- -->

``` r
# ggsave(here("figs", "revision_figs", "figsupp2_bulk_density.png"), width=200, height=160, units="mm", dpi=400)
```

``` r
bd_min_max <- bd_depth %>%
  dplyr::summarize(min = min(depth_wt_bd, na.rm=TRUE),
            max = max(depth_wt_bd, na.rm=TRUE))
```

``` r
bd_soc <- esm %>%
  left_join(dplyr::select(bd_depth, -depth), by=c("dsp_pedon_id", "project", "depth_std" = "layer", "label"))

treatment_bd_cv <- bd_soc %>%
  filter(method=="fd", project!="Illinois", project!="UConn") %>%
  group_by(project, label, depth_std) %>%
  dplyr::summarize(mean_bd = mean(depth_wt_bd),
            sd_bd = sd(depth_wt_bd)) %>%
  mutate(cv_bd_treatment = sd_bd/mean_bd)
```

    ## `summarise()` has grouped output by 'project', 'label'. You can override using
    ## the `.groups` argument.

``` r
overall_treatment_cv <- treatment_bd_cv %>%
  ungroup() %>%
  group_by(depth_std, label) %>%
  dplyr::summarize(treatment_mean_cv = mean(cv_bd_treatment, na.rm=TRUE))
```

    ## `summarise()` has grouped output by 'depth_std'. You can override using the
    ## `.groups` argument.

``` r
project_bd_cv <- bd_soc %>%
  filter(method=="fd", project!="Illinois", project!="UConn") %>%
  group_by(project, depth_std) %>%
  dplyr::summarize(mean_bd = mean(depth_wt_bd),
            sd_bd = sd(depth_wt_bd)) %>%
  mutate(cv_bd_project = sd_bd/mean_bd)
```

    ## `summarise()` has grouped output by 'project'. You can override using the
    ## `.groups` argument.

``` r
overall_project_cv <- project_bd_cv %>%
  ungroup() %>%
  group_by(depth_std) %>%
  dplyr::summarize(project_mean_cv = mean(cv_bd_project, na.rm=TRUE))

bd_cv_table <- overall_treatment_cv %>%
  left_join(overall_project_cv, by="depth_std") %>%
  mutate(across(treatment_mean_cv:project_mean_cv, ~round(.x, 2))) %>%
  arrange(depth_std, label)

# write_csv(bd_cv_table, here("figs", "revision_figs", "table1_bd_cv.csv"))
```

## Range of reference masses for different soils

``` r
mean_ref_mass <- esm %>%
  group_by(depth_std) %>%
  dplyr::summarize(mean_mass = mean(soil_mass),
            mean_mass_cum = mean(soil_mass_cum))

cust_labeller <- function(x) paste0("ESM Depth: ", x)

ggplot(esm %>% filter(method!="fd"), aes(x=project, y=soil_mass, fill=method)) +
  geom_boxplot() +
  facet_wrap(~esm_depth, scales="free", labeller=as_labeller(cust_labeller)) +
  scale_fill_paletteer_d("nationalparkcolors::Arches",
                         name="Calculation method") + 
  labs(x="DSP4SH project", y="Reference soil mass (Mg/ha)") +
  scale_x_discrete(labels=project_labels_esm) +
  theme_classic() +
    theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="none")
```

![](esm_manuscript_revision_files/figure-gfm/fig%202%20range%20of%20reference%20masses%20in%20each%20ESM%20depth-1.png)<!-- -->

``` r
# ggsave(here("figs", "revision_figs", "fig2_soil_mass.png"), width=180, height=160, units="mm", dpi=400)
# ggsave(here("figs", "revision_figs", "fig2.pdf"), width=180, height=160, units="mm", dpi=400)
```

``` r
ref_mass_table <- esm %>%
  group_by(depth_std) %>%
  dplyr::summarize(min_mass = min(soil_mass),
            max_mass = max(soil_mass))
```

# How does ESM choice influence total calculated SOC stocks?

## Effect of stock calculation method on cumulative SOC stocks

What is the mean SOC stock calculated via each method?

What is the difference between the smallest and largest SOC stocks
calculated for a particular pedon using these different methods?

``` r
soc_method_summary <- esm %>%
  group_by(depth_std, method_longest) %>%
  dplyr::summarize(mean = round(mean(soc_cum),1))
```

    ## `summarise()` has grouped output by 'depth_std'. You can override using the
    ## `.groups` argument.

``` r
soc_method_minmaxdiff <- esm %>%
  filter(depth_std=="60") %>%
  group_by(project, label, dsp_pedon_id, depth_std) %>%
  dplyr::summarize(min_soc_cum = min(soc_cum),
            max_soc_cum = max(soc_cum),
            diff_soc_cum = max_soc_cum - min_soc_cum,
            pct_diff_soc_cum = ((max_soc_cum - min_soc_cum) / min_soc_cum)*100)
```

    ## `summarise()` has grouped output by 'project', 'label', 'dsp_pedon_id'. You can
    ## override using the `.groups` argument.

``` r
mean_diff_project <- soc_method_minmaxdiff %>%
  ungroup() %>%
  group_by(project) %>%
  dplyr::summarize(mean_diff_soc_cum = round(mean(diff_soc_cum), 2),
            min_diff_soc_cum = round(min(diff_soc_cum), 2),
            max_diff_soc_cum = round(max(diff_soc_cum), 2))

flextable(mean_diff_project)
```

<img src="esm_manuscript_revision_files/figure-gfm/soc stock methods summary-1.png" width="1199" />

``` r
mean_diff <- soc_method_minmaxdiff %>%
  ungroup() %>%
  dplyr::summarize(mean_diff_soc_cum = round(mean(diff_soc_cum), 2),
            min_diff_soc_cum = round(min(diff_soc_cum), 2),
            max_diff_soc_cum = round(max(diff_soc_cum), 2),
            max_pct_diff_soc_cum = round(max(pct_diff_soc_cum), 1),
            mean_pct_diff_soc_cum = round(mean(pct_diff_soc_cum), 1))

flextable(mean_diff)
```

<img src="esm_manuscript_revision_files/figure-gfm/soc stock methods summary-2.png" width="1557" />

Run linear mixed effects model:

``` r
soc_method_lmer <- esm %>%
  group_by(depth_std) %>%
  nest() %>%
  mutate(lmer = purrr::map(data, ~lmer(soc_cum ~ method_longest + (1|project/label), data = .x)),
  drop1 = purrr::map(data, .f = ~{
    drop1(lmer(soc_cum ~ method_longest + (1|project/label), REML=FALSE, data = .x))
    }),
  letters = purrr::map(lmer, .f = ~{
    glht <- glht(.x, linfct = mcp(method_longest = "Tukey"))
    cld <- cld(glht)
    tidy(cld)
    }),
  tidy_lmer = purrr::map(lmer, broom.mixed::tidy),
  tidy_drop = purrr::map(drop1, broom.mixed::tidy),
  predict = purrr::map(lmer, ~ggpredict(.x, terms = c("method_longest")))) 
```

    ## Warning: There were 4 warnings in `mutate()`.
    ## The first warning was:
    ## ℹ In argument: `tidy_drop = purrr::map(drop1, broom.mixed::tidy)`.
    ## ℹ In group 1: `depth_std = 5`.
    ## Caused by warning:
    ## ! The column names NumDF and DenDF in ANOVA output were not recognized or
    ## transformed.
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 3 remaining warnings.

``` r
# Extract table of F-values and p-values for each depth
method_drop_table <- soc_method_lmer %>%
  select(depth_std, tidy_drop) %>%
  unnest(cols=c(tidy_drop), names_sep="_") %>%
  mutate(sig = case_when(tidy_drop_p.value < 0.05 ~ "significant",
                         tidy_drop_p.value > 0.05 ~ "not_significant")) %>%
  select(depth_std, tidy_drop_statistic, tidy_drop_p.value, sig) %>%
  distinct() %>%
  mutate(across(where(is.numeric), ~round(.x, 2)),
         tidy_drop_p.value = ifelse(tidy_drop_p.value<0.001, "<0.001", tidy_drop_p.value))

flextable(method_drop_table)
```

<img src="esm_manuscript_revision_files/figure-gfm/lmer calculation method-1.png" width="866" />

``` r
# Extract significantly different management groups for each method
method_cld_table <- soc_method_lmer %>%
  select(depth_std, letters) %>%
  unnest(cols=c(letters))

flextable(method_cld_table)
```

<img src="esm_manuscript_revision_files/figure-gfm/lmer calculation method-2.png" width="653" />

These plots show that (unsurprisingly), using different calculation
methods results in significantly different SOC stocks, especially in
surface soil horizons. Using ESM min methods results in the smallest
stocks (of course), while fixed depth tends to be more comparable to ESM
mean estimates.

``` r
fig2_letters_position <- esm %>%
  group_by(depth_std, method_longest) %>%
  dplyr::summarize(max_soc_cum = max(soc_cum)) %>%
  left_join(method_cld_table, by=c("depth_std", "method_longest"))
```

    ## `summarise()` has grouped output by 'depth_std'. You can override using the
    ## `.groups` argument.

``` r
ggplot(esm, aes(x=method_longest, y=soc_cum, fill=method_longest)) +
  geom_point(aes(color=method_longest), position = "jitter", alpha=0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_text(data=fig2_letters_position, aes(y=max_soc_cum + 30, x=method_longest, label=letters), size=3) +
  labs(x="Calculation method", y="Cumulative SOC stock (Mg/ha)") +
  facet_wrap(~depth_std, labeller=labeller(depth_std=cum_depth_labels)) +
  scale_x_discrete(labels=method_labels) +
  scale_color_paletteer_d("nationalparkcolors::Arches",
                          name="Calculation method") + 
  scale_fill_paletteer_d("nationalparkcolors::Arches",
                          name="Calculation method") + 
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="none")
```

![](esm_manuscript_revision_files/figure-gfm/fig%203%20plot%20actual%20cumulative%20SOC%20stocks%20at%20each%20depth-1.png)<!-- -->

``` r
# ggsave(here("figs", "revision_figs", "fig3_soc_cum_all.png"), width=160, height=120, units="mm", dpi=400)
# ggsave(here("figs", "revision_figs", "fig3.pdf"), width=160, height=120, units="mm", dpi=400)
```

# How does calculation method influence management sensitivity?

## Calculate means for all treatments, depths, and calculation methods for results text

``` r
soc_cum_summary <- esm %>%
  group_by(method_longest, esm_depth, label) %>%
  dplyr::summarize(mean = round(mean(soc_cum), 1)) %>%
  arrange(esm_depth, label, method_longest)
```

    ## `summarise()` has grouped output by 'method_longest', 'esm_depth'. You can
    ## override using the `.groups` argument.

``` r
soc_cum_summary2 <- soc_cum_summary %>%
  group_by(esm_depth, label) %>%
  dplyr::summarize(min = min(mean),
            max = max(mean))
```

    ## `summarise()` has grouped output by 'esm_depth'. You can override using the
    ## `.groups` argument.

## Management sensitivity of cumulative SOC stocks

We can test the overall sensitivity to management by running a mixed
effects linear model that accounts for depth and project:

``` r
soc_mgmt_lmer <- esm %>%
  group_by(method_longest, depth_std) %>%
  nest() %>%
  mutate(lmer = purrr::map(data, ~lmer(soc_cum ~ label + (1|project), data = .x)),
  drop1 = purrr::map(data, .f = ~{
    drop1(lmer(soc_cum ~ label + (1|project), REML=FALSE, data = .x))
    }),
  letters = purrr::map(lmer, .f = ~{
    glht <- glht(.x, linfct = mcp(label = "Tukey"))
    cld <- cld(glht)
    tidy(cld)
    }),
  tidy_lmer = purrr::map(lmer, broom.mixed::tidy),
  tidy_drop = purrr::map(drop1, broom.mixed::tidy),
  predict = purrr::map(lmer, ~ggpredict(.x, terms = c("label")))) 
```

    ## Warning: There were 20 warnings in `mutate()`.
    ## The first warning was:
    ## ℹ In argument: `tidy_drop = purrr::map(drop1, broom.mixed::tidy)`.
    ## ℹ In group 1: `method_longest = esm2_min_project` `depth_std = 5`.
    ## Caused by warning:
    ## ! The column names NumDF and DenDF in ANOVA output were not recognized or
    ## transformed.
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 19 remaining warnings.

``` r
# Extract table of F-values and p-values for each method
mgmt_drop_table <- soc_mgmt_lmer %>%
  select(method_longest, depth_std, tidy_drop) %>%
  unnest(cols=c(tidy_drop), names_sep="_") %>%
    select(method_longest, tidy_drop_statistic, tidy_drop_p.value) %>%
  mutate(tidy_drop_p.value = ifelse(tidy_drop_p.value < 0.05,  "<0.001", as.character(tidy_drop_p.value)),
         tidy_drop_statistic = round(tidy_drop_statistic, 2),
         method_longest=factor(method_longest, levels=method_factor_order)) %>%
  arrange(method_longest, depth_std) %>%
  relocate(method_longest, .before=depth_std)
```

    ## Adding missing grouping variables: `depth_std`

``` r
# Extract significantly different management groups for each method
mgmt_cld_table <- soc_mgmt_lmer %>%
  select(method_longest, depth_std, letters) %>%
  unnest(cols=c(letters))

tablesupp2 <- mgmt_drop_table %>%
  left_join(mgmt_cld_table, by=c("method_longest", "depth_std")) %>%
  mutate(label = factor(label, levels=c("BAU", "SHM", "Ref")),
         method_longest=factor(method_longest, levels=method_factor_order)) %>%
  arrange(method_longest, depth_std, label)

flextable(tablesupp2)
```

<img src="esm_manuscript_revision_files/figure-gfm/lmer management-1.png" width="1345" />

``` r
# write_csv(tablesupp2, here("figs", "revision_figs", "tablesupp3_mgmt.csv"))
```

When accounting for differences in depth and project, all calculation
methods result in SOC stocks that are highest in Ref condition, and
lower in both SHM and BAU conditions.

Plot actual differences due to management:

``` r
mgmt_letters_df <- esm %>%
  group_by(method_longest, depth_std) %>%
  dplyr::summarize(max_soc_cum = max(soc_cum)) %>%
  left_join(mgmt_cld_table, by=c("depth_std", "method_longest")) %>%
  mutate(label = factor(label, levels=c("BAU", "SHM", "Ref")))
```

    ## `summarise()` has grouped output by 'method_longest'. You can override using
    ## the `.groups` argument.

``` r
ggplot(esm, aes(x=method_longest, y=soc_cum, fill=label)) +
  geom_point(aes(color=label), pch = 21, position=position_jitterdodge(), alpha=0.8) +
  geom_boxplot(outlier.shape=NA) +
  scale_x_discrete(labels=method_labels) +
  facet_wrap(~depth_std, scales="free", labeller=labeller(depth_std=cum_depth_labels)) +
  geom_text(data=mgmt_letters_df, 
            aes(y=max_soc_cum + 25, label=letters, group=label), 
            color="black", position=position_dodge(width=0.8), size=3) +
  labs(x="Calculation method", y="Cumulative SOC (Mg/ha)") +
  scale_fill_manual(values=c("#FED789FF","#72874EFF","#476F84FF"),
                    breaks=c("BAU", "SHM", "Ref"), 
                    name="Management") +
  scale_color_manual(values=c("#FED789FF","#72874EFF","#476F84FF"),
                     breaks=c("BAU", "SHM", "Ref"), 
                     name="Management") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
```

![](esm_manuscript_revision_files/figure-gfm/fig%204%20plot%20actual%20cumulative%20SOC%20stocks%20under%20different%20management-1.png)<!-- -->

``` r
# ggsave(here("figs", "revision_figs", "fig4_soc_mgmt.png"), width=180, height=160, units="mm", dpi=400)
# ggsave(here("figs", "revision_figs", "fig4.pdf"), width=180, height=160, units="mm", dpi=400)
```

# Calculation method influence on SOC sequestration

Calculate the mean SOC in the BAU treatment for each project, and then
calculate the difference between each Ref/SHM pedon and the BAU mean.

This isn’t entirely different from assessing the management sensitivity
of each calculation method, as we did above. However, in this analysis
we are concerned with the magnitude of the difference in SOC stocks
between both Ref and SHM treatments compared to BAU.

``` r
esm_bau <- esm %>%
  filter(label=="BAU") %>%
  group_by(project, method_long, method_longest, method, ref_stat, ref_data, depth_std) %>%
  dplyr::summarize(mean_bau_soil_mass = mean(soil_mass),
            mean_bau_soc = mean(soc),
            mean_bau_soil_mass_cum = mean(soil_mass_cum),
            mean_bau_soc_cum = mean(soc_cum))
```

    ## `summarise()` has grouped output by 'project', 'method_long', 'method_longest',
    ## 'method', 'ref_stat', 'ref_data'. You can override using the `.groups`
    ## argument.

``` r
esm_diff <- esm %>%
  filter(!label=="BAU") %>%
  left_join(esm_bau, by=c("project", "method_long", "method_longest", "method", "ref_stat", "ref_data", "depth_std")) %>%
  mutate(diff_soil_mass = soil_mass - mean_bau_soil_mass,
         diff_soc = soc - mean_bau_soc,
         diff_soil_mass_cum = soil_mass_cum - mean_bau_soil_mass_cum,
         diff_soc_cum = soc_cum - mean_bau_soc_cum)

# Calculate mean sequestration estimates for each treatment, depth, and method
esm_diff_summary <- esm_diff %>%
  group_by(label, depth_std, method_longest) %>%
  dplyr::summarize(mean_diff = round(mean(diff_soc_cum), 1),
            n = n())
```

    ## `summarise()` has grouped output by 'label', 'depth_std'. You can override
    ## using the `.groups` argument.

``` r
flextable(esm_diff_summary)
```

<img src="esm_manuscript_revision_files/figure-gfm/calculate difference between BAU and Ref/SHM SOC stocks-1.png" width="965" />

## Cumulative SOC sequestration with different calculation methods

``` r
diff_cum_method_lmer <- esm_diff %>%
  group_by(depth_std, label) %>%
  nest() %>%
  mutate(lmer = purrr::map(data, ~lmer(diff_soc_cum ~ method_longest + (1|project), data = .x)),
  drop1 = purrr::map(data, .f = ~{
    drop1(lmer(diff_soc_cum ~ method_longest + (1|project), REML=FALSE, data = .x))
    }),
  letters = purrr::map(lmer, .f = ~{
    glht <- glht(.x, linfct = mcp(method_longest = "Tukey"))
    cld <- cld(glht)
    tidy(cld)
    }),
  tidy_lmer = purrr::map(lmer, broom.mixed::tidy),
  tidy_drop = purrr::map(drop1, broom.mixed::tidy),
  predict = purrr::map(lmer, ~ggpredict(.x, terms = c("method_longest")))) 
```

    ## Warning: There were 8 warnings in `mutate()`.
    ## The first warning was:
    ## ℹ In argument: `tidy_drop = purrr::map(drop1, broom.mixed::tidy)`.
    ## ℹ In group 1: `depth_std = 5` `label = SHM`.
    ## Caused by warning:
    ## ! The column names NumDF and DenDF in ANOVA output were not recognized or
    ## transformed.
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 7 remaining warnings.

``` r
# Extract table of F-values and p-values for each depth
diff_cum_method_drop_table <- diff_cum_method_lmer %>%
  unnest(cols=c(tidy_drop), names_sep="_") %>%
  mutate(sig = case_when(tidy_drop_p.value < 0.05 ~ "significant",
                         tidy_drop_p.value > 0.05 ~ "not_significant")) %>%
  select(depth_std, label, tidy_drop_statistic, tidy_drop_p.value, sig) %>%
  distinct() %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

flextable(diff_cum_method_drop_table)
```

<img src="esm_manuscript_revision_files/figure-gfm/lmer for influence of calculation method on differences in cumulative soc stocks-1.png" width="1076" />

``` r
# Extract significantly different management groups for each method
diff_cum_method_cld_table <- diff_cum_method_lmer %>%
  select(depth_std, label, letters) %>%
  unnest(cols=c(letters))

flextable(diff_cum_method_cld_table)
```

<img src="esm_manuscript_revision_files/figure-gfm/lmer for influence of calculation method on differences in cumulative soc stocks-2.png" width="803" />

``` r
diff_cum_letters_position <- esm_diff %>%
  group_by(method_longest, depth_std, label) %>%
  dplyr::summarize(max_diff_soc_cum = max(diff_soc_cum)) %>%
  left_join(diff_cum_method_cld_table, by=c("depth_std", "label", "method_longest"))
```

    ## `summarise()` has grouped output by 'method_longest', 'depth_std'. You can
    ## override using the `.groups` argument.

``` r
ggplot(esm_diff, aes(x=method_longest, y=diff_soc_cum, fill=method_longest)) +
  geom_point(aes(color=method_longest), pch=21, position="jitter", alpha=0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_abline(intercept=0, slope=0, linetype="dashed") +
  geom_text(data=diff_cum_letters_position, aes(y=max_diff_soc_cum + 25, label=letters), size=3) +
  labs(x="Calculation method", y="SOC stock difference vs BAU (Mg/ha)") +
  facet_grid(depth_std~label, scales="free", labeller=labeller(depth_std=cum_depth_labels)) +
  scale_x_discrete(labels=method_labels) +
  scale_fill_paletteer_d("nationalparkcolors::Arches",
                         name="Calculation method") + 
  scale_color_paletteer_d("nationalparkcolors::Arches",
                         name="Calculation method") + 
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="none")
```

![](esm_manuscript_revision_files/figure-gfm/fig%205%20plot%20actual%20differences%20in%20cumulative%20SOC%20stocks%20with%20different%20methods-1.png)<!-- -->

``` r
# ggsave(here("figs", "revision_figs", "fig5_soc_diff.png"), width=180, height=170, units="mm", dpi=400)
# ggsave(here("figs", "revision_figs", "fig5.pdf"), width=180, height=170, units="mm", dpi=400)
```

# Check that ESM calculations worked and make sense

## Plot mean cumulative SOC vs depth, and cumulative SOC vs cumulative soil mass

These plots are just to check that the ESM calculations worked and
nothing looks too crazy - for example, wild fluctuations in SOC content
with depth, impossible numbers, etc.

Cumulative SOC vs depth:

``` r
purrr::map(.x = projects, 
           .f = ~{
             esm %>%
               filter(project==.x) %>%
               ggplot(aes(x=depth_std, y=soc_cum, color=method_longest)) +
               stat_summary(geom="line", fun="mean", linewidth=1.3) +
               coord_flip() +
               xlim(c(60,0)) +
               facet_grid(~label, scales="free") +
               labs(x="ESM Depth (cm)", y="Cumulative SOC (Mg/ha)", 
                    title=glue::glue("Cumulative SOC with depth - ", .x)) +
               scale_color_paletteer_d("nationalparkcolors::Arches",
                                       name="Calculation method", labels=method_labels) + 
               theme_classic()
           })
```

    ## [[1]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20depth-1.png)<!-- -->

    ## 
    ## [[2]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20depth-2.png)<!-- -->

    ## 
    ## [[3]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20depth-3.png)<!-- -->

    ## 
    ## [[4]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20depth-4.png)<!-- -->

    ## 
    ## [[5]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20depth-5.png)<!-- -->

    ## 
    ## [[6]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20depth-6.png)<!-- -->

    ## 
    ## [[7]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20depth-7.png)<!-- -->

    ## 
    ## [[8]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20depth-8.png)<!-- -->

    ## 
    ## [[9]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20depth-9.png)<!-- -->

Cumulative SOC vs cumulative soil mass:

``` r
purrr::map(.x = projects, 
           .f = ~{
             esm %>%
               filter(project==.x) %>%
               ggplot(aes(x=soil_mass_cum, y=soc_cum, color=label, group=dsp_pedon_id)) +
               geom_smooth(stat="identity") +
               coord_flip() +
               scale_x_reverse() +
               facet_wrap(~method_longest, labeller=labeller(method_longest=method_labels2)) +
               labs(x="Cumulative soil mass (Mg/ha)", y="Cumulative SOC (Mg/ha)", 
                    title=glue::glue("Cumulative SOC vs cumulative soil mass - ", .x)) +
               scale_color_manual(values=c("#FED789FF","#72874EFF","#476F84FF"),
                                  breaks=c("BAU", "SHM", "Ref"), 
                                  name="Management") +
               theme_classic()
           })
```

    ## [[1]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20cumulative%20soil%20mass-1.png)<!-- -->

    ## 
    ## [[2]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20cumulative%20soil%20mass-2.png)<!-- -->

    ## 
    ## [[3]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20cumulative%20soil%20mass-3.png)<!-- -->

    ## 
    ## [[4]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20cumulative%20soil%20mass-4.png)<!-- -->

    ## 
    ## [[5]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20cumulative%20soil%20mass-5.png)<!-- -->

    ## 
    ## [[6]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20cumulative%20soil%20mass-6.png)<!-- -->

    ## 
    ## [[7]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20cumulative%20soil%20mass-7.png)<!-- -->

    ## 
    ## [[8]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20cumulative%20soil%20mass-8.png)<!-- -->

    ## 
    ## [[9]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soc%20vs%20cumulative%20soil%20mass-9.png)<!-- -->

These aren’t the easiest to interpret, but I think that everything looks
good here. It doesn’t look like SOC estimates take a crazy left turn at
any depth for any of the methods, and they are all in relative agreement
with each other. The plot of cumulative mass vs cumulative SOC shows
that using the minimum reference mass for ESM results in a lower
cumulative soil mass overall, but that’s no surprise. Fixed depth
results look pretty identical to ESM mean.

## How does ESM choice influence soil mass in each depth increment?

This is another way to check that nothing is getting wonky with the ESM
calculations - making sure that the soil mass in each depth increment
makes sense for each project.

Plot soil mass in each depth increment:

``` r
# Soil mass in each depth increment
purrr::map(.x = projects, 
           .f = ~{
             esm %>%
               filter(project==.x) %>%
               ggplot(aes(x=esm_depth, y=soil_mass, fill=method_longest)) +
               stat_summary(geom="bar", fun=mean, position="dodge2", color="black") +
               stat_summary(geom="errorbar", fun.data = mean_cl_normal, color="black",
                            width=0.4, position=position_dodge(width=0.9)) +
               coord_flip() +
               scale_x_discrete(limits=rev) +
               facet_wrap(~label) +
               labs(x="ESM Depth", y="Depth Increment Soil Mass (Mg/ha)", 
                    title=glue::glue("Depth increment soil mass by calculation method - ", .x)) +
               scale_fill_paletteer_d("nationalparkcolors::Arches",
                                      name="Calculation method",
                                      labels=method_labels) + 
               theme_classic()
           })
```

    ## [[1]]

![](esm_manuscript_revision_files/figure-gfm/soil%20mass%20in%20each%20depth%20increment-1.png)<!-- -->

    ## 
    ## [[2]]

![](esm_manuscript_revision_files/figure-gfm/soil%20mass%20in%20each%20depth%20increment-2.png)<!-- -->

    ## 
    ## [[3]]

![](esm_manuscript_revision_files/figure-gfm/soil%20mass%20in%20each%20depth%20increment-3.png)<!-- -->

    ## 
    ## [[4]]

![](esm_manuscript_revision_files/figure-gfm/soil%20mass%20in%20each%20depth%20increment-4.png)<!-- -->

    ## 
    ## [[5]]

![](esm_manuscript_revision_files/figure-gfm/soil%20mass%20in%20each%20depth%20increment-5.png)<!-- -->

    ## 
    ## [[6]]

![](esm_manuscript_revision_files/figure-gfm/soil%20mass%20in%20each%20depth%20increment-6.png)<!-- -->

    ## 
    ## [[7]]

![](esm_manuscript_revision_files/figure-gfm/soil%20mass%20in%20each%20depth%20increment-7.png)<!-- -->

    ## 
    ## [[8]]

![](esm_manuscript_revision_files/figure-gfm/soil%20mass%20in%20each%20depth%20increment-8.png)<!-- -->

    ## 
    ## [[9]]

![](esm_manuscript_revision_files/figure-gfm/soil%20mass%20in%20each%20depth%20increment-9.png)<!-- -->

Plot cumulative soil mass in each depth increment:

``` r
# Cumulative soil mass in each depth increment
purrr::map(.x = projects, 
           .f = ~{
             esm %>%
               filter(project==.x) %>%
               ggplot(aes(x=factor(depth_std), y=soil_mass_cum, fill=method_longest)) +
               stat_summary(geom="bar", fun=mean, position="dodge2", color="black") +
               stat_summary(geom="errorbar", fun.data=mean_cl_normal, color="black",
                            width=0.4, position=position_dodge(width = 0.9)) +
               coord_flip() +
               scale_x_discrete(limits=rev) +
               facet_wrap(~label) +
               labs(x="ESM Depth (cm)", y="Cumulative Soil Mass (Mg/ha)", 
                    title=glue::glue("Cumulative soil mass by calculation method - ", .x)) +
               scale_fill_paletteer_d("nationalparkcolors::Arches",
                                      name="Calculation method",
                                      labels=method_labels) + 
               theme_classic()
           })
```

    ## [[1]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soil%20mass%20in%20each%20depth%20increment-1.png)<!-- -->

    ## 
    ## [[2]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soil%20mass%20in%20each%20depth%20increment-2.png)<!-- -->

    ## 
    ## [[3]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soil%20mass%20in%20each%20depth%20increment-3.png)<!-- -->

    ## 
    ## [[4]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soil%20mass%20in%20each%20depth%20increment-4.png)<!-- -->

    ## 
    ## [[5]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soil%20mass%20in%20each%20depth%20increment-5.png)<!-- -->

    ## 
    ## [[6]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soil%20mass%20in%20each%20depth%20increment-6.png)<!-- -->

    ## 
    ## [[7]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soil%20mass%20in%20each%20depth%20increment-7.png)<!-- -->

    ## 
    ## [[8]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soil%20mass%20in%20each%20depth%20increment-8.png)<!-- -->

    ## 
    ## [[9]]

![](esm_manuscript_revision_files/figure-gfm/cumulative%20soil%20mass%20in%20each%20depth%20increment-9.png)<!-- -->

These plots are a check to see how using ESM methods influence the mass
of soil in a given depth increment. They show that, unsurprisingly,
fixed depth calculations result in the greatest soil mass in each depth
increment. Fixed depth soil masses are most comparable to ESM(mean)
masses.
