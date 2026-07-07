# Helper functions for ESM calculations :)

# Function to calculate aggregated masses using hybrid bulk density in each depth increment ----
# Returns minimum and average values
# Required input: dataframe of soil horizon data, vector with desired bottom depths of each depth increment
soil_mass_aggregate_hybrid <- function(input, depth){
  # Promote to SPC and dice into 1-cm increments
  spc <- input
  aqp::depths(spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
  
  # dice into 1-cm intervals
  dice <- aqp::dice(spc, fm=0:99 ~ bd_hybrid)
  
  dice_mass <- horizons(dice) %>%
    mutate(mass = (hrzdep_b - hrzdep_t) * bd_hybrid * 100)
  
  #Initialize output vector
  out <- vector("list", length(depth))
  
  min_mean_max <- list(
    min = ~min(.x, na.rm=TRUE), 
    mean = ~mean(.x, na.rm=TRUE),
    max = ~max(.x, na.rm=TRUE)
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
    dplyr::summarize(across(mass_agg, min_mean_max))
  
}

# Function to calculate SOC stocks AND SOIL MASSES with depth increments that mirror ESM increments ----
soc_stock_fd_hybrid <- function(input, depth){
  spc <- input
  aqp::depths(spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
  
  # dice into 1-cm intervals
  dice <- aqp::dice(spc, fm=0:99 ~ bd_hybrid + soc_fill)
  
  soc <- horizons(dice) %>%
    mutate(hrzdepth = hrzdep_b - hrzdep_t) %>%
    mutate(soc_stock_hrz = soc_fill * bd_hybrid * hrzdepth,
           soil_mass_hrz = bd_hybrid * hrzdepth * 100)
  
  soc_out <- vector("list", length(depth))
  
  for (i in seq_along(depth)) {
    if (i == 1) {
      soc_out[[i]] <- soc %>%
        filter(hrzdep_t < depth[[i]]) %>%
        group_by(dsp_pedon_id) %>%
        dplyr::summarize(soc.fd = sum(soc_stock_hrz)) %>%
        mutate(topdepth.fd = 0,
               depth_cat = depth[[i]])
      
    } else {
      soc_out[[i]] <- soc %>%
        filter(hrzdep_t >=depth[[i-1]] & hrzdep_t < depth[[i]]) %>%
        group_by(dsp_pedon_id) %>%
        dplyr::summarize(soc.fd = sum(soc_stock_hrz)) %>%
        mutate(topdepth.fd = -depth[[i-1]],
               depth_cat = depth[[i]])
    }
  }
  
  mass_out <- vector("list", length(depth))
  
  for (i in seq_along(depth)) {
    if (i == 1) {
      mass_out[[i]] <- soc %>%
        filter(hrzdep_t < depth[[i]]) %>%
        group_by(dsp_pedon_id) %>%
        dplyr::summarize(soil_mass.fd = sum(soil_mass_hrz)) %>%
        mutate(topdepth.fd = 0,
               depth_cat = depth[[i]])
      
    } else {
      mass_out[[i]] <- soc %>%
        filter(hrzdep_t >=depth[[i-1]] & hrzdep_t < depth[[i]]) %>%
        group_by(dsp_pedon_id) %>%
        dplyr::summarize(soil_mass.fd = sum(soil_mass_hrz)) %>%
        mutate(topdepth.fd = -depth[[i-1]],
               depth_cat = depth[[i]])
    }
  }
  
  soc_df <- dplyr::bind_rows(soc_out) %>%
    group_by(dsp_pedon_id) %>%
    mutate(layer=seq_along(depth_cat)) %>%
    mutate(depth.fd = -depth_cat) %>%
    select(-depth_cat) %>%
    unite("sample_id", c("dsp_pedon_id", "layer"), sep="-", remove=FALSE) %>%
    ungroup()
  
  mass_df <- dplyr::bind_rows(mass_out) %>%
    group_by(dsp_pedon_id) %>%
    mutate(layer=seq_along(depth_cat)) %>%
    select(-depth_cat) %>%
    ungroup()
  
  soc_mass_df <- soc_df %>%
    left_join(mass_df, by=c("dsp_pedon_id", "topdepth.fd", "layer")) %>%
    group_by(dsp_pedon_id) %>%
    arrange(dsp_pedon_id, desc(depth.fd)) %>%
    mutate(soil_mass_cum.fd = cumsum(soil_mass.fd),
           soc_cum.fd = cumsum(soc.fd)) %>%
    relocate(soc.fd, .before=soil_mass.fd)
  
  soc_mass_df
  
}
