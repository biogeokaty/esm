# Helper functions for ESM calculations :)

# Function to calculate estimated horizon boundaries for generalized horizons ----
# Requirements: SPC object, vector with generalized horizon names, vector with generalized horizon patterns

genhz_depths <- function(spc, names, patt){
  # Add generalized horizon labels, remove unused levels, and make factor
  spc$genhz <- generalize.hz(spc$hzdesg, names, patt) 
  spc$genhz[spc$genhz == "not-used"] <- NA
  spc$genhz <- factor(spc$genhz)
  
  # keep track of generalized horizon names for later
  hz_names <- levels(spc$genhz)
  
  # define vector for dicing
  slice_vect <- seq(from = 0, to = 99, by = 1)
  
  # dice into 1-cm intervals
  dice <- aqp::dice(spc, slice_vect ~ genhz)
  dice$genhz <- factor(dice$genhz, levels = hz_names)
  
  # Logistic Proportional-Odds Ordinal Regression Model
  dd <- rms::datadist(horizons(dice))
  options(datadist = "dd")
  l.genhz <- orm(genhz ~ rcs(hrzdep_t), data = horizons(dice), x = TRUE, y = TRUE)
  
  # predict along same depths: columns are the class-wise probability fitted.ind --> return all probability estimates
  predict <- data.frame(predict(l.genhz, data.frame(hrzdep_t = slice_vect), type = "fitted.ind"))
  
  # re-name, rms model output give funky names
  names(predict) <- hz_names
 
   # add depths
  predict$top <- slice_vect
  
  # maximum likelihood depths
  predict_ml <- get.ml.hz(predict, o.names = hz_names)
  
  # return output
  predict_ml
}

# Function to calculate aggregated masses in each depth increment and return max, min, and average values ----
soil_mass_aggregate <- function(input, depth){
  out <- vector("list", length(depth))
  min_max_mean <- list(
    min = ~min(.x, na.rm=TRUE), 
    max = ~max(.x, na.rm=TRUE),
    mean = ~mean(.x, na.rm=TRUE)
  )
  
  for (i in seq_along(depth)) {
    if (i == 1) {
      out[[i]] <- input %>%
        filter(hrzdep_t < depth[[i]]) %>%
        group_by(dsp_pedon_id) %>%
        dplyr::summarize(mass_agg = sum(mass)) %>%
        mutate(depth_cat = depth[[i]])
      
    } else {
      out[[i]] <- input %>%
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

# Function to calculate SOC stocks with depth increments that mirror ESM increments ----
soc_stock_fd <- function(input, depth){
  
  soc <- aqp::horizons(input) %>%
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
