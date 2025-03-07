# Helper functions for ESM calculations :)

# Function to promote SPC object, assign generalized horizon names, and calculate depth increments ----
# Input: dataframe object, vector with generalized horizon names, vector with generalized horizon patterns, vector of desired dicing depth increments

genhz_depths <- function(dataframe, names, patt){

  # Promote to SPC
  spc <- dataframe
  aqp::depths(spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
  aqp::hzdesgnname(spc) <- 'hzdesg'
  
  # Add generalized horizon labels, remove unused levels, and make factor
  spc$genhz <- aqp::generalize.hz(spc$hzdesg, names, patt) 
  spc$genhz[spc$genhz == "not-used"] <- NA
  spc$genhz <- factor(spc$genhz)
  
  # keep track of generalized horizon names for later
  hz_names <- levels(spc$genhz)
  
  # dice into 1-cm intervals
  dice <- aqp::dice(spc, seq(from = 0, to = 99, by = 1) ~ genhz)
  dice$genhz <- factor(dice$genhz, levels = hz_names)
  
  # Logistic Proportional-Odds Ordinal Regression Model
  horizons <- aqp::horizons(dice)
  dd <- rms::datadist(horizons)
  options(datadist = "dd")
  l.genhz <- orm(genhz ~ rcs(hrzdep_t), data = horizons, x = TRUE, y = TRUE)
  
  # predict along same depths: columns are the class-wise probability fitted.ind --> return all probability estimates
  predict <- data.frame(predict(l.genhz, data.frame(hrzdep_t = seq(from = 0, to = 99, by = 1)), type = "fitted.ind"))
  
  # re-name, rms model output give funky names
  names(predict) <- hz_names
  
  # add depths
  predict$top <- seq(from = 0, to = 99, by = 1)
  
  # maximum likelihood depths
  predict_ml <- get.ml.hz(predict, o.names = hz_names)
  
  # return output
  predict_ml
}

# Function to calculate aggregated masses in each depth increment and return max, min, and average values ----
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

# Function to calculate SOC stocks with depth increments that mirror ESM increments ----
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
        mutate(topdepth_fd = 0,
               depth_cat = depth[[i]])
      
    } else {
      out[[i]] <- soc %>%
        filter(hrzdep_t >=depth[[i-1]] & hrzdep_t < depth[[i]]) %>%
        group_by(dsp_pedon_id) %>%
        dplyr::summarize(soc_fd = sum(soc_stock_hrz)) %>%
        mutate(topdepth_fd = -depth[[i-1]],
               depth_cat = depth[[i]])
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
