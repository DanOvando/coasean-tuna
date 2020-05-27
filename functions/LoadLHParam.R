
### Function that loads life history and stock assessment data by species. 

LoadLHParam <- function(species_names,
                        data_year,
                        fleet_names,
                        results_dir){
  
  # Path to data
  dat_dir <- here::here("data", paste0(data_year, "_SA"), paste0(species_names, "/"))
  pr_dat_dir <- here::here("data", paste0("price/"))
  
  ## Place to save output .csv files and data plots
  save_dir <- paste0(results_dir, data_year, "/", "inputs/")
  
  if (dir.exists(save_dir) == F) {
    dir.create(save_dir, recursive = T)
  }
  
  # Load general life history characteristics
  dat <- read_csv(here::here("data", paste0(data_year, "_SA"), paste0("life_history_inputs_",data_year,".csv"))) %>%
    dplyr::filter(species %in% species_names) %>%
    arrange(species)

  # Make life history list
  lh <- NULL
  lh$species_name <- dat$species
  lh$max_age <- dat$max_age
  ages <- seq(1:lh$max_age)
  
  # Natural mortality
  m <- read_delim(paste0(dat_dir, species_names, "_m_", data_year, ".txt"), 
                  delim = '\t')$y
  
  if(length(m) == 1){
    lh$m <- rep(m, length = length(ages))
  }else if(length(m) == length(ages)){
    lh$m <- m 
  }else{
    lh$m <- rep(m, length = length(ages))
    warning("Number of supplied natural mortality values does not match number of age classes")
  }
  
  # Natural survival
  lh$s <- exp(-lh$m)
  
  # Length at age
  lh$length_at_age <- read_delim(paste0(dat_dir, species_names, "_length_at_age_", data_year, ".txt"),
                      delim = '\t')$y
  
  if(length(lh$length_at_age) != length(ages)){
    warning("Number of supplied lengths-at-age values does not match number of age classes")
  }
  
  # Weight at age
  lh$wa <- dat$wa # parameter a for L vs W
  lh$wb <- dat$wb # parameter b for L vs W
  lh$weight_form <- as.character(dat$weight_form)
  lh$weight_sd <- dat$sd_weight
  lh$weight_sd_limits <- c(dat$sd_weight_low, dat$sd_weight_high)
  
  if(lh$weight_form == "Exponential"){
    
    weights <- (lh$wa*(lh$length_at_age^lh$wb))
    
  }else if(lh$weight_form == "Linear"){
    
    weights <- (lh$wa+(l_age*lh$wb))
    
  }else{
    warning("Not a valid weight form")
  }
  
  lh$weight_at_age <- weights * 0.001 # convert kg to tons
  
  # Maturity at age
  lh$maturity_at_age <- read_delim(paste0(dat_dir, species_names, "_mat_", data_year, ".txt"),
                                   delim = '\t')$y
  
  if(length(lh$maturity_at_age) != length(ages)){
    warning("Number of supplied maturity values does not match number of age classes")
  }
  
  # Life history plots
  
  lh_df <- tibble(species = species_names,
                  age = ages,
                  m = lh$m,
                  length_at_age = lh$length_at_age,
                  weight_at_age = lh$weight_at_age,
                  maturity_at_age = lh$maturity_at_age)
  
  m_plot <- ggplot(lh_df, aes(x = age, y = m)) +
    geom_line(color = "black") +
    theme_bw()
  
  length_at_age_plot <- ggplot(lh_df, aes(x = age, y = length_at_age)) +
    geom_line(color = "red") +
    theme_bw()
  
  weight_at_age_plot <- ggplot(lh_df, aes(x = age, y = weight_at_age)) +
    geom_line(color = "red") +
    annotate("text", x = min(lh_df$age), y = max(lh_df$weight_at_age)-0.001, label = paste('W = ', lh$wa, '*L^(', lh$wb, ')', sep = ''), hjust = 0)+
    theme_bw()
  
  maturity_at_age_plot <- ggplot(lh_df, aes(x = age, y = maturity_at_age)) +
    geom_line(color = "black") +
    theme_bw()
  
  lh_plots <- length_at_age_plot + weight_at_age_plot + m_plot + maturity_at_age_plot 
  
  ggsave(paste0(save_dir, "life_history_", species_names, ".png"), lh_plots, dpi = 200)
  
  # BH Recruitment
  lh$sex_ratio <- 0.5
  lh$BH_steepness <- dat$BH_steepness
  
  # Estimated virgin population values
  lh$B0 <- NA 
  lh$R0 <- NA
  lh$SSB0 <- NA
  lh$SSB0_R0 <- NA
  
  # Estimated stock reference points
  lh$Bmsy<- NA
  lh$Fmsy <- NA
  lh$Fmsy_tot <- NA
  lh$SSBmsy <- NA
  lh$MSY <- NA
  
  # Comparison reference points from the stock assessment 
  ref_points <- read_csv(here::here("data", paste0(data_year, "_SA"), paste0("ref_points_", data_year, ".csv"))) %>%
    dplyr::filter(species %in% species_names) %>%
    arrange(species)
  
  lh$ref_points_known <- ref_points
  lh$B0_known <- ref_points$B0
  lh$SSB0_known <- ref_points$SSB0
  lh$SSBmsy_known <- ref_points$SSBmsy
  lh$Fmsy_tot_known <- ref_points$Fmsy
  lh$MSY_known <- ref_points$msy
  
  ### Depletion (SB/SBf=0)
  
  depletion_df <- read_delim(paste0(dat_dir, species_names, "_dep_", data_year, ".txt"), 
                  delim = '\t') %>%
    rename(depletion = y)
  
  depletion <- as.matrix(depletion_df %>% dplyr::select(-x))
  rownames(depletion) <- round(depletion_df$x, digits = 0)
  
  lh$depletion_avg <- depletion
  
  ### Fishing Mortality 
  
  juv_f_df <- read_delim(paste0(dat_dir, species_names, "_f_juv_", data_year, ".txt"), 
                         delim = '\t') %>%
    rename(juv = y)
  
  adult_f_df <- read_delim(paste0(dat_dir, species_names, "_f_adult_", data_year, ".txt"), 
                         delim = '\t') %>%
    rename(adult = y)
  
  f_df <- juv_f_df %>%
    bind_cols(adult_f_df)
  
  fishing_mort <- as.matrix(f_df %>% dplyr::select(-x, -x1))
  rownames(fishing_mort) <- round(adult_f_df$x, digits = 0)
  
  lh$f_avg <- fishing_mort
  
  ### Recruitment (input from SA is in millions of fish)
  
  rec_df <- read_delim(paste0(dat_dir, species_names, "_recruitment_", data_year, ".txt"), 
                       delim = '\t') %>%
    mutate(year = round(x, digits = 0)) %>%
    mutate(recruits = y*1000000)
  
  recruits <- as.matrix(rec_df %>% dplyr::select(-year, -y, -x))
  rownames(recruits) <- rec_df$year
  
  lh$recruitment_total <- recruits
  
  rec_start <- min(as.numeric(rownames(recruits)))
  rec_end <- max(as.numeric(rownames(recruits)))
  
  recruitment_series <- approx(rec_df$year, rec_df$recruits/4, method = "constant", seq(rec_start, (rec_end + 0.75), by = 0.25), rule = 2)$y %>%
    as.matrix()
  rownames(recruitment_series) <- seq(rec_start, (rec_end + 0.75), by = 0.25)
  
  lh$recruitment_qtr <- recruitment_series
  
  ### SSB (input from SA is in 1000's of mt)

  ssb_df <- read_delim(paste0(dat_dir, species_names, "_ssb_", data_year, ".txt"), 
                       delim = '\t') %>%
    mutate(ssb = y * 1000)
  
  ssb <- as.matrix(ssb_df %>% dplyr::select(-x, -y))
  rownames(ssb) <- round(ssb_df$x, digits = 0)
  
  lh$ssb_avg <- ssb 
  
  # ### Prices -----------------

  price_dat <- read_csv(paste0(pr_dat_dir, "price_dat.csv"))
  
  price_blank <- tibble(fleet = rep(fleet_names, each = length(ages)),
                        age_quarters = rep(ages, times = length(fleet_names)),
                        species = species_names) %>%
    mutate(gear_type = case_when(grepl("PS", fleet) ~ "PS",
                                 grepl("PL", fleet) ~ "PL",
                                 grepl("HL", fleet) ~ "HL",
                                 grepl("LL", fleet) ~ "LL",
                                 TRUE ~ "OTH"))
  
  prices_fleet <- price_blank %>%
    left_join(price_dat %>% dplyr::filter(species == species_names), by = c("gear_type", "age_quarters", "species"))
  
  price_mat <- as.matrix(prices_fleet %>%
                           dplyr::select(fleet, age_quarters, price) %>%
                            group_by(fleet) %>%
                            spread(age_quarters, price) %>%
                            ungroup() %>%
                            arrange(fleet) %>%
                            dplyr::select(-fleet))
  rownames(price_mat) <- fleet_names
  
  lh$prices <- price_mat

# Return
  return(lh)
}
