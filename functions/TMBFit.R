###
### Script that uses TMB to fit R0 and reference points for an age-structured population with multiple fleets. 
### TMB is estimating qs over time for each fleet in order to explain catches and relative depletion. 
###

TMBFit <- function(lh,
                   catch,
                   selec,
                   data_year,
                   fleet_names,
                   effort,
                   results_dir,
                   historic_weight = 0.5,
                   version = "msy",
                   use_historic_recruits = FALSE){
  
### General values ----------------------
  
  species_names <- lh$species_name

  if(data_year == 2014){
    years <- c(1952, 2012)
  }else if(data_year == 2017){
    years <- c(1952, 2015)
  }

  save_dir <- paste0(results_dir, data_year, "/", "bio_fits/")
  
  if (dir.exists(save_dir) == F) {
    dir.create(save_dir, recursive = T)
  }
  
### Ensure data is in correct formats -----------------------
  
  ### Depletion data (year)
  model_years <- rownames(lh$ssb_avg)
  
  ssb_series <- subset(lh$ssb_avg, subset = rownames(lh$ssb_avg) %in% model_years)[,1]
  # subset(lh$ssb_avg, subset = rownames(lh$ssb_avg) %in% model_years)[,1]
  
  model_years <- as.numeric(names(ssb_series))
  
  model_quarters <- seq(first(model_years), last(model_years) + 0.75, by = 0.25) # in year.quarter format 
  
  ## depletion
  depletion_series <- subset(lh$depletion_avg, subset = rownames(lh$depletion_avg) %in% model_years)[,1]
  
  ### Recruitment data (quarter)
  recruitment_series <- subset(lh$recruitment_qtr, subset = rownames(lh$recruitment_total) %in% model_quarters)[,1]
  
  if (use_historic_recruits == FALSE){
    recruitment_series <- rep(0, length(recruitment_series))
    
  }

  recruitment_quarters <- sum(is.na(recruitment_series) == F)
  ### Selectivity data (fleets x age_classes)
  selectivity_dat <- selec

  ### Catch data (model_quarters x fleets)
  catch_dat <- catch
  
  catch_years <- floor(as.numeric(rownames(catch_dat))) %in% model_years
  
  catch_dat <- catch_dat[catch_years, ] 
  
  ### Effort data (model_quarters x fleets)
  effort_dat <- effort
  
  effort_years <- floor(as.numeric(rownames(effort_dat))) %in% model_years
  
  effort_dat <- effort_dat[effort_years, ] 
  
  any_catch <- catch_dat > 0
  
  effort_dat <- effort_dat * any_catch
  
### Assign TMB input data ------------------------------------------------------------------------------------------

 weights <- c(0.9, 0.1) # SSB, catch

project_recruit_data <- rep(mean(recruitment_series), length(recruitment_series))

 data <- list(max_age = lh$max_age,
              quarters = length(model_quarters),
              years = length(model_years),
              num_fleets = length(fleet_names),
              steepness = lh$BH_steepness,
              data_start_quarter = 0,
              m = lh$m,
              s = lh$s,
              weight_at_age = lh$weight_at_age,
              maturity_at_age = lh$maturity_at_age,
              selectivities = selectivity_dat %>% t(),
              catches = catch_dat,
              efforts = effort_dat,
              recruit_data = recruitment_series,
              project_recruit_data = project_recruit_data,
              recruit_quarters = recruitment_quarters,
              depletion = depletion_series, # length years
              dep_weight = weights,
              msy = lh$MSY_known,
              fmult = seq(0, 5, length.out = 50),
              ssb_msy_to_ssb_0 = lh$SSBmsy_known/lh$SSB0_known)
 

 params <- list(log_r0 = log(6e9),
                 q_est = matrix(log(1e-6), length(unique(rownames(catch_dat))) *
                                  length(unique(colnames(catch_dat))),
                                nrow = data$num_fleets,
                                ncol = data$quarters) %>% t(),
                log_burn_f = log(mean(lh$m)/3), # f during burn in period (pre-catches),
                log_dep_sigma = log(0.05)) # ln(sd_dev) for depletion
                # log_catch_sigma = log(1),
                # log_recruit_sigma = log(1)) # ln(sd_dev) for catch
  
  q_length <- length(unique(rownames(catch_dat))) * length(unique(colnames(catch_dat)))


#### ----------------------------------------------------------------------------------------------------------
### Use the new TMB model which DOES try to calculate and fit to MSY/other reference points    
### DAN - troubleshooting this model is above my capabilities at the momement. You rewrote this script in order to be able to 
### return multiple outputs, and it would take me more time than I have available at the moment to figure this out. 
  
  # browser()
  # wtf <- params
  compile(here::here("scripts", paste0("tuna_model_",version,".cpp")), flags = "-std=c++11") 
  dyn.load(dynlib(here::here("scripts", paste0("tuna_model_",version))))
  model <- MakeADFun(data, params, DLL=paste0("tuna_model_",version))
  lower = c(log(1e-20), rep(log(1e-20), q_length), rep(log(1e-20), length(model$par) - q_length - 1)) # lower bound for all model parameters by position
  upper = c(log(1e20), rep(log(1e-1), q_length), rep(log(1), length(model$par) - q_length - 1)) # upper bound for all model parameters by position
  
  ### Fit TMB model------------------------------------------------------------------------------------------
  
  seed <- round(runif(1, min = 100, max = 999))
  set.seed(seed)
  fit <- nlminb(model$par, model$fn, model$gr, control = list(iter.max=5000, eval.max=5000),
                lower = lower, upper = upper)
  
  # for (i in 1:3) # do it 3 more times just to make sure
  #   fit <- nlminb(model$env$last.par.best , model$fn, model$gr, control = list(iter.max=1000, eval.max = 1000),
  #                 lower = lower, upper = upper)
  
  fit_ss <- fit$objective
  
  fit_report <- model$report()
  
  fit_report$msy_location <-  fit_report$msy_location + 1

### Reformat outputs ------------------------------------------------------------------------------------------
  
  pop <- model$report()$n_at_a
  rownames(pop) <- model_quarters
  colnames(pop) <- seq(1, lh$max_age, 1)
  
  biomass <- model$report()$b_at_a
  rownames(biomass) <- model_quarters
  colnames(biomass) <- seq(1, lh$max_age, 1)
  
  SSB_age <- model$report()$ssb_at_a
  rownames(SSB_age) <- model_quarters
  
  SSB <- rowSums(SSB_age) %>% as.matrix()
  rownames(SSB) <- model_quarters
  
  recruits <- model$report()$recruits %>% as.matrix()
  rownames(recruits) <- model_quarters

  f_age <- model$report()$f_age
  rownames(f_age) <- model_quarters
  colnames(f_age) <- seq(1, lh$max_age, 1)
  
  mean_f_age <- rowMeans(f_age) %>% as.matrix()
  rownames(mean_f_age) <- model_quarters
  
  f_total <- rowSums(f_age) %>% as.matrix()
  rownames(f_total) <- model_quarters
  
  f_fleet_apical <- model$report()$f_fleet_apical
  colnames(f_fleet_apical) <- fleet_names
  rownames(f_fleet_apical) <- model_quarters
  
  q_fleet <- model$report()$q_at_t
  colnames(q_fleet) <- fleet_names
  rownames(q_fleet) <- model_quarters
  
  catches_age <- model$report()$c_at_a
  rownames(catches_age) <- model_quarters

  catch_total <- rowSums(catches_age) %>% as.matrix()
  rownames(catch_total) <- model_quarters

  catches_fleet <- model$report()$c_at_fleet
  colnames(catches_fleet) <- fleet_names
  rownames(catches_fleet) <- model_quarters
  
  # Initial conditions
  lh$pop <- pop
  lh$SSB_age <- SSB_age
  
  lh$f_fleet_recent <- fit_report$f_fleet_recent
  lh$fmult_catches <- fit_report$fmult_catches
  lh$fmult_ssb <- fit_report$fmult_ssb
  
  lh$catches_fleet <- catches_fleet
  lh$catch_total <- catch_total
  lh$catch_ss <- model$report()$c_ss
  lh$dep_ss <- model$report()$dep_ss
  
  lh$B0 <- model$report()$b0
  lh$N0 <- model$report()$n0
  lh$R0 <- model$report()$r0
  lh$SSB0 <- model$report()$total_ssb0
  lh$SSB0_R0 <- lh$SSB0/lh$R0
  lh$alpha <- model$report()$alpha
  lh$beta <- model$report()$beta
  lh$depletion_avg_hat <- model$report()$dep_avg
  
  # Reference points
  lh$fmult <- fit_report$fmult[fit_report$msy_location]
  lh$Bmsy <- fit_report$ssb_msy # FIX THIS
  lh$Nmsy <- fit_report$ssb_msy # FIX THIS
  lh$SSBmsy <- fit_report$ssb_msy # FIX THIS
  lh$MSY <- fit_report$msy_hat
  
  fmsy <- fit_report$fmsy
  names(fmsy) <- fleet_names
  
  lh$Fmsy <- sum(fmsy) # FIX THIS - should be sum of fs by age NOT FLEET
  lh$Fmsy_fleet <- fmsy 
  lh$Fmsy_meanage <- mean(fit_report$fmsy) # FIX THIS - should be average f by age NOT FLEET
  
  # Yield multipler plot 
  yield_mult_plot <- tibble(multiplier = fit_report$fmult,
                                yield = fit_report$fmult_catches) %>% 
    ggplot() + 
    aes(multiplier, yield)+
    geom_point() + 
    geom_vline(aes(xintercept = fit_report$fmult[fit_report$msy_location]), linetype = 2, color = 'red')
  
  ggsave(paste0(save_dir, "yield_mult_plot_", species_names, ".png"), yield_mult_plot, dpi = 200)
  
### Make Kobe plot  ------------------------------------------------------------------------------------------
  
  
  SSB_time_series <- SSB %>%
    as.tibble() %>%
    mutate(year_quarter = as.numeric(rownames(SSB))) %>%
    rename(SSB = V1) %>%
    mutate(year = floor(year_quarter)) %>%
    group_by(year) %>%
    summarize(SSB = mean(SSB)) %>%
    mutate(SSB_SSBmsy = SSB/lh$SSBmsy)
  
  lh$f_fleet_apical <-  f_fleet_apical
  
  f_time_series <- tibble(year_quarter = as.numeric(rownames(f_fleet_apical)),
                          f = rowSums(f_fleet_apical)) %>% 
    mutate(year = floor(year_quarter)) %>%
    group_by(year) %>%
    summarize(f = mean(f)) %>%
    mutate(F_Fmsy = f/lh$Fmsy)
  
  
    # f_fleet_apical %>%
    # rowSums()
    # as.tibble() %>%
    # mutate(year_quarter = as.numeric(rownames(f_age))) %>%
    # gather(fleet, f_age, -year_quarter) %>%
    # group_by(year_quarter) %>%
    # summarize(f_age_mean = mean(f_age)) %>%
    # ungroup() %>%
    # mutate(year = floor(year_quarter)) %>%
    # group_by(year) %>%
    # summarize(f = mean(f_age_mean)) %>%
    # mutate(F_Fmsy = f/lh$Fmsy)
  
  time_series <- SSB_time_series %>%
    left_join(f_time_series, by = "year") %>%
    mutate(species = species_names)
  
  lh$ref_points <- time_series
  # xxmax <- max(2, ceiling(1.4*max(time_series$SSB_SSBmsy)))
  # yymax <- max(2, ceiling(1.4*max(time_series$F_Fmsy)))
  # 
  # kobe_plot <- ggplot(time_series, aes(x = SSB_SSBmsy, y = F_Fmsy, fill = year))+
  #   theme_bw()+
  #   theme(plot.margin = unit(c(0.5, 0.5 , 0.5, 0.5), "cm"),
  #         axis.ticks = element_line(color = "black", size = 1),
  #         axis.ticks.length = unit(0.5, "cm"),
  #         axis.text = element_text(size = 14),
  #         axis.title = element_text(size = 16),
  #         legend.background = element_blank(),
  #         legend.key = element_blank())+
  #   geom_rect(xmin = 1, xmax = xxmax, ymin = 0, ymax = 1, fill = 'green')+
  #   geom_rect(xmin = 0, xmax = 1, ymin = 1, ymax = yymax, fill = 'red')+
  #   geom_rect(xmin = 1, xmax = xxmax, ymin = 1, ymax = yymax, fill = 'orange')+
  #   geom_rect(xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = 'yellow')+
  #   geom_hline(yintercept = 1, size = 1)+
  #   geom_vline(xintercept = 1, size = 1)+
  #   scale_y_continuous(expand = c(0.015,0), limits = c(0, yymax), breaks = c(seq(from = 0, to = yymax, by = 1)))+
  #   scale_x_continuous(expand = c(0.015,0), limits = c(0, xxmax), breaks = c(seq(from = 0, to = xxmax, by = 1)))+
  #   geom_point(shape = 21, size = 7) +
  #   geom_path(color = "black")+
  #   scale_fill_distiller(palette = "Purples", direction = 1, limits = c(min(time_series$year), max(time_series$year)))+
  #   labs(x= "SSB/SSBmsy", y="F/Fmsy", title = paste(species_name, " Kobe Plot SSB/SSBmsy vs F/Fmsy (total)", sep = ""))+
  #   geom_segment(size = 1, x = 0, xend=xxmax, y = -0.025, yend = -0.025)+
  #   geom_segment(size = 1, x = -0.045, xend= -0.045, y = 0, yend=yymax)
  # 
  # ggsave(paste0(save_dir, "kobe_plot_", species_names, ".png"), kobe_plot, dpi = 200)
 
  
  return(lh)
  
}