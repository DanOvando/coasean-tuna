# -------------------------------------------------------------------------
# Run Coaseian Tuna
# A wrapper for running the the Coaseian tuna analysis
# Also generates figures 3 and 4 for the paper and the figures for the supplementary materials
# -------------------------------------------------------------------------

# Load packages
set.seed(123)

# Analysis
library(ranger)
library(caret)
library(gbm)
library(broom)
library(pbapply)
library(doMC)
library(lme4)
library(TMB)
library(recipes)
library(parsnip)
library(patchwork)
library(rsample)


# Plotting
library(ggplot2)
library(scales) # Scales
library(sf) # Spatial work
library(viridis) # Color palettes
library(ggpubr) # Multiple plot arragements

#library(patchwork)

# General
library(here)
library(tidyverse)
# library(extrafont)
# extrafont::loadfonts()
# library(hrbrthemes)


functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions

# Parameters --------------------------------------------------------------

run_name <-  'mre-submission'

run_description <-
  "PNA only version of paper submission for MRE"

assessment_year <-  2014

recent_year <- 2011

fit_models <- FALSE

run_biological_fits <- FALSE

fit_bargain <-  FALSE

run_in_parallel <-  TRUE

num_cores <- 8

num_chains <- 1

smooth_things <- F

whale_sharks <- T

remove_oth <-  T

pna_only <- T

iterations <- 6000

species_of_interest = c('skipjack', 'bigeye') #species to deal with, ignore yellowfin basically

cut_year <- 1950 #cutoff year for data to include

reg_year <- 1995 #cutoff year for data to include in regression

npv_factor <-
  1e9 # express NPV in terms of what, default to billions

skipjack_price = 1782 # skipjack  Price per unit weight(?)

# skipjack_price = 100 # skipjack  Price per unit weight(?)


bigeye_price = 9370 # targeted bigeye price per unit weight(?)

bigeye_bycatch_price = 1782 # bycatch bigeye price per unit weight

manual_bet_fad_effect = 0.688 # from summer month comparison in explore_tuna_data

manual_skj_fad_effect = 0.048 # from summer month comparison in explore_tuna_data

disc_rate <- 0.1

manual_fad_effect <-  F

use_stan <- F

use_glmer <- T

use_hurdle <- T

form = 'hurdle'

litraining_set = 'all'

run_dir = paste('results/', run_name, '/', sep = '')

dep_var = 'log_cpue'

drop_countries = c('New Zealand', 'South Korea', 'Japan', 'China')

ind_vars = c('fad',
             'factor_month',
             'factor_country',
             # 'factor_year',
             'factor_area',
             'sst',
             'sst2')

# tree_vars <-  c('fad',
#                 'factor_area',
#                 'factor_country',
#                 'sst',
#                 'month',
#                 'year')

tree_vars <-  c('fad',
                'lon',
                'lat',
                'sst',
                'month',
                'year',
                'factor_country')



reg_fmla <-
  'norm_dep_var ~ (1 + fad | factor_area:factor_month) + factor_country  + sst + sst2'

# reg_fmla <- 'norm_dep_var ~ (1 + fad | factor_area:factor_month) + factor_country + factor_year  + sst + sst2'

# reg_fmla <- NA

logit_fmla <- NA

gam_smooth = c('days')

drop_months = 999

if (dir.exists(run_dir) == F) {
  dir.create(run_dir, recursive = T)
  
}

write(run_description, file = paste(run_dir, 'DESCRIPTION.txt', sep = ''))

# Set Plot Themes --------------------------------------------------------------

title_size <-
  12 # Plot titles (axis, legends, plot, etc) - 14 for paper, 20 for ppts
text_size <-
  10 # All other text (axis text, legend text, etc.) - 12 for paper, 18 for ppts
point_size <- 5

simple_theme <- theme_bw() +
  theme(
    text = element_text(size = text_size),
    title = element_text(size = title_size),
    legend.position = "bottom"
  )

theme_set(simple_theme)

# Export width for figures (in) - this applies to all figures
fig_width = 7.5

# Expore height for standard figures (in) - this only applies to misc. figures, all publication figures have height specified as a proportion of width individually
fig_height = 4.5

# Load  catch Data ---------------------------------------------------------------
# lookup tables of gear, fleet, and species for various indices

gear_lookup = read.csv(file = here('data', 'gear_lookup.csv'),
                       stringsAsFactors = F) %>%
  mutate(gear_name = tolower(gear_name))

fleet_lookup = read.csv(file = here('data', 'fleet_lookup.csv'),
                        stringsAsFactors = F) %>%
  mutate(fleet_name = tolower(fleet_name))

species_lookup = read.csv(file = here('data', 'species_lookup.csv'),
                          stringsAsFactors = F) %>%
  mutate(species_name = tolower(species_name),
         species = tolower(species))

# Load in raw spatial data  ------------------------------------------------

small_num <-  1e-6

min_catch <-  -1

sst_data <-
  read.csv(here('data','SST' ,'wcpfc_sst_00_15.csv'), stringsAsFactors = F) %>%
  dplyr::select(-contains('X')) %>%
  rename(
    lon = LON5,
    lat = LAT5,
    year = YY,
    month = MM,
    sst = avg_sst
  ) %>%
  as_tibble() %>%
  unique()

eez_data = read.csv(here('data', 'EEZ summaries', 'PS_5_y_EEZ_detailed.csv'),
                    stringsAsFactors = F) %>%
  as_tibble() %>%
  rename(
    lon = LON5_D,
    lat = LAT5_D,
    eez = EEZ,
    country = Country,
    pna_country = PNA_country
  ) %>%
  dplyr::select(lon, lat, eez, country, pna_country) %>%
  ungroup() %>%
  unique()


# tidy purse seine catches --------------

ps_dbf <-
  foreign::read.dbf(here("data", "catch_effort", "purse_seine_4", "PURSE_SEINE.DBF")) %>%
  set_names(tolower(colnames(.))) %>%
  as_tibble() %>%
  rename(year = yy, month = mm)

purse_data <- ps_dbf %>%
  gather(catch_type, catch, contains('_c_')) %>%
  select(-contains('sets')) %>%
  mutate(
    species = str_split(catch_type, '_c_', simplify = T)[, 1],
    set_type = str_split(catch_type, '_c_', simplify = T)[, 2]
  ) %>% {
    if (whale_sharks == T) {
      mutate(., fad = set_type %>% as.factor())
    } else {
      mutate(., fad = str_detect(set_type, '(fad)|(log)'))
      
    }
    
  } %>%
  group_by(year, month, lat5, lon5, species, fad) %>%
  summarise(catch = sum(catch), days = unique(days))

# tidy set allocation data ----------------

set_data <- ps_dbf %>%
  gather(set_type, sets, contains('sets')) %>%
  select(-contains('_c_')) %>%
  mutate(set_type = str_split(set_type, '_', simplify = T)[, 2]) %>% {
    if (whale_sharks == T) {
      mutate(., fad = set_type %>% as.factor())
    } else {
      mutate(., fad = str_detect(set_type, '(fad)|(log)'))
      
    }
    
  } %>%
  group_by(year, month, lat5, lon5, fad) %>%
  summarise(sets = sum(sets)) %>%
  group_by(year, month, lat5, lon5) %>%
  mutate(total_sets = sum(sets),
         prop_sets = sets / total_sets)

# add set data to catch data -------------

purse_data <- purse_data %>%
  left_join(set_data, by = c('year', 'month', 'lat5', 'lon5', 'fad')) %>%
  ungroup() %>%
  mutate(days = days * prop_sets) %>%
  filter(days > 0 & is.na(days) == F)

# add purse seine covariates ------------


purse_data <- purse_data %>%
  mutate(
    north_south = str_sub(lat5, start = -1, end = -1),
    #convert to decimal degrees
    east_west = str_sub(lon5, start = -1, end = -1),
    ns_index = ifelse(north_south == 'N', 1, -1),
    ew_index = ifelse(east_west == 'E', 1, -1),
    lat = ns_index * as.numeric(gsub('.{1}$', '', lat5)),
    lon = ew_index * as.numeric(gsub('.{1}$', '', lon5))
  ) %>%
  left_join(sst_data, by = c('lat', 'lon', 'month', 'year')) %>%
  left_join(eez_data, by = c('lat', 'lon'))  # add in eez data

sst_model <-
  lm(sst ~ lon + lat + lon * lat + month %>% as.factor(), data = sst_data)

pred_sst <- predict(sst_model,
                    newdata = purse_data)

purse_data$sst[is.na(purse_data$sst)] <-
  pred_sst[is.na(purse_data$sst)]

purse_data$sst2 <- purse_data$sst ^ 2


# process purse_seine covariates --------

# new region definitions that need to be worked in:
# region 1: 20 - 50 lat and 120 - 170 lon
# region 2: 20 - 50 lat and 170 - 210 lon
# region 3: 0 - 20 lat and 140 - 170 lon, or -5 - 0 lat and 155 - 170 lon, or -10 t0 - 5 lat and 160 - 170 lon
# region 4: -10 to 20 lat and 170 - 210 lon
# region 5: -15 to -10 lat and 140 - 170 lon, or -20 to -15 lat and 150 to 170 lon, or -45 to -20 lat and 140 - 170 lon
# region 6: -40 to -10 lat and 170 to 210 lon
# region 7: -10 to 20 lat and 110 to 140 lon
# region 8: -10 to 0 lat and 140 to 155 lon, or -10 to -5 lat and 155 to 160 lon
# region 9: -20 to -15 lat and 140 to 150 lon
# area = ifelse(lat >= 0,
#               ifelse(lon >= 140, 2, 1),
#               ifelse(lon < 140, 3, 4))


purse_data <- purse_data %>%
  mutate(
    #Calculate summary statistics
    area = ifelse(
      lon < 0,
      ifelse(lat > 20, 2,
             ifelse(lat %in% -9:20, 4, 6)),
      
      ifelse(
        lat > 20,
        ifelse(lon > 170, 2, 1),
        ifelse(
          lat %in% 1:20,
          ifelse(lon <= 140, 7,
                 ifelse(lon > 170, 4, 3)),
          
          ifelse(
            lat %in% -4:0,
            ifelse(lon <= 140, 7,
                   ifelse(
                     lon %in% 141:155, 8,
                     ifelse(lon %in% 156:170, 3, 4)
                   )),
            
            ifelse(
              lat %in% -9:-5,
              ifelse(lon <= 140, 7,
                     ifelse(
                       lon %in% 141:160, 8,
                       ifelse(lon %in% 161:170, 3, 4)
                     )),
              
              ifelse(
                lat %in% -14:-10,
                ifelse(lon <= 170, 5, 6),
                
                ifelse(
                  lat %in% -19:-15,
                  ifelse(lon <= 150, 9,
                         ifelse(lon <= 170, 5, 6)),
                  
                  ifelse(lat <= -20,
                         ifelse(lon <= 170, 5, 6), "bad")
                )
              )
            )
          )
        )
      )
    ),
    
    log_catch = log(catch + small_num),
    log_cpue = log((catch / days + small_num)),
    cpue = catch / days,
    power_days = days ^ 1.5,
    # power_days2 = days ^ 1.5,
    factor_country = as.factor(country),
    factor_month = as.factor(month),
    factor_area = as.factor(area),
    any_caught = catch > min(catch),
    factor_year = as.factor(year),
    size_bin = cut(log_catch, 4)
  ) %>%
  mutate(entry = paste(year, month, lat, lon, sep = '-'))

ggplot(purse_data, aes(x = lon, y = lat, color = factor_area)) +
  geom_point()

purse_data$factor_country <-
  as.character(purse_data$factor_country)

purse_data$factor_country[is.na(purse_data$factor_country)] <-
  'unknown'

purse_data$factor_country <-  factor(purse_data$factor_country)

# Filter data -------------------------------------------------------------

purse_data <- purse_data %>%
  filter(!(country %in% drop_countries))

# Run fad removal regression ------------------------------------------

bet_dat = purse_data %>%
  filter(species == 'bet' &
           year >= cut_year & sets > 0 & is.na(sets) == F) %>%
  group_by(factor_country) %>%
  mutate(num_obs = sum(cpue > 0 & is.na(cpue) == F)) %>%
  filter(num_obs > 10) %>%
  ungroup()

skj_dat = purse_data %>%
  filter(species == 'skj' &
           year >= cut_year  & sets > 0 & is.na(sets) == F) %>%
  group_by(factor_country) %>%
  mutate(num_obs = sum(cpue > 0 & is.na(cpue) == F)) %>%
  filter(num_obs > 10) %>%
  ungroup()

partition_tuna <- function(i, dat, test_set) {
  if (i == 'all') {
    training <-  dat
    
    test <- dat
    
  } else {
    training <- slice(dat,-test_set[i])
    
    test <- slice(dat, test_set[i])
    
  }
  
  out = list(training = training, test = test)
  
  return(out)
}

ind_vars_2 <- ind_vars[ind_vars != 'factor_month']

center_scale <- function(x, xname, omit_names = '') {
  if (is.numeric(x) & !all(unique(x) %in% c(1, 0)) &
      !xname %in% omit_names) {
    x <- (x - mean(x, na.rm = T)) / (2 * sd(x, na.rm = T))
    
  }
  
  return(x)
  
}

# center and scale data
# purse_data <-  purse_data %>%
#   map2_df(colnames(.), ~center_scale(.x,.y, omit_names = .y[!(.y %in% ind_vars)])) %>%
#   mutate(sst2 = sst ^2)





if (fit_models == T) {
  
  bet_dat <- tibble::rownames_to_column(bet_dat)
  
  skj_dat <- tibble::rownames_to_column(skj_dat)
  
  # tree_vars <- c(tree_vars,"rowname")
  
  bet_model_comp <- fit_forests(
    dat = bet_dat,
    reg_year = reg_year,
    tree_vars = tree_vars,
    ind_vars = ind_vars,
    glmer_fmla = reg_fmla,
    dep_var = dep_var,
    form = form,
    use_hurdle = use_hurdle,
    iterations = iterations,
    run_in_parallel = run_in_parallel,
    num_cores = num_cores
  )
  
  
  skj_model_comp <- fit_forests(
    dat = skj_dat,
    reg_year = reg_year,
    tree_vars = tree_vars,
    ind_vars = ind_vars,
    glmer_fmla = reg_fmla,
    dep_var = dep_var,
    form = form,
    iterations = iterations,
    use_hurdle = use_hurdle,
    run_in_parallel = run_in_parallel,
    num_cores = num_cores
    
  )
  
  bet_model_comp_future <- fit_forests(
    dat = bet_dat,
    reg_year = reg_year,
    tree_vars = tree_vars,
    ind_vars = ind_vars,
    glmer_fmla = reg_fmla,
    dep_var = dep_var,
    form = form,
    use_hurdle = use_hurdle,
    iterations = iterations,
    run_in_parallel = run_in_parallel,
    num_cores = num_cores,
    holdout_mode = 2
  )
  
  
  skj_model_comp_future <- fit_forests(
    dat = skj_dat,
    reg_year = reg_year,
    tree_vars = tree_vars,
    ind_vars = ind_vars,
    glmer_fmla = reg_fmla,
    dep_var = dep_var,
    form = form,
    iterations = iterations,
    use_hurdle = use_hurdle,
    run_in_parallel = run_in_parallel,
    num_cores = num_cores,
    holdout_mode = 2
  )
  
  write_rds(skj_model_comp_future, file.path(run_dir, "skj_model_comp_future.rds"))
  
  write_rds(bet_model_comp_future, file.path(run_dir, "bet_model_comp_future.rds"))
  
  
  save(file = paste0(run_dir, 'coaseian_trees.Rdata'),
       bet_model_comp,
       skj_model_comp)
} else{
  load(paste0(run_dir, 'coaseian_trees.Rdata'))
  
  skj_model_comp_future <- read_rds(file.path(run_dir, "skj_model_comp_future.rds"))
  
  bet_model_comp_future <-  read_rds(file.path(run_dir, "bet_model_comp_future.rds"))
  
  
  
  # load(paste(run_dir, 'coaseian_fad_bargain_results.Rdata', sep = ''))
  
  functions <- list.files(here::here("functions"))
  
  walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions
  #save_plots <- TRUE
  
  # load_functions(func_dir = 'functions')
  
  # run_dir = paste('results/', run_name, '/', sep = '')
  
  
  
} #close if fit model

message('done with model fitting')


### Appendix Figures - Model fits ### -----------------------------------------------------

### Figure A1 - Predicted vs. Observed CPUE for both BET and SKJ
# FAD closure time frame(s)
closed_frame <- bet_model_comp$model_comp %>%
  mutate(numyear = as.numeric(year_month)) %>%
  select(year, month, year_month) %>%
  filter(month == 7, year > 2008) %>%
  mutate(date = lubridate::ymd(paste(year, month, "1", sep = '/')))

consistent_bet_places <- bet_dat %>% 
  filter(year < 2009,fad != "una") %>% 
  select(month,lat, lon) %>% 
  unique() %>% 
  mutate(mll = paste(month,lat,lon, sep = '-'))

# BET

bet_fad_effect <- bet_model_comp$forest_predictions %>%
  mutate(cpue = (catch / days),
         date = lubridate::ymd(paste(year, month, "1", sep = '/')),
         mll = paste(month,lat,lon, sep = '-')) %>%
  filter(year > 2008,
         mll %in% consistent_bet_places$mll) %>%
  gather("cpue_type","cpue", cpue, cpue_hat) %>% 
  group_by(date, cpue_type) %>%
  summarise(mean_cpue = mean(cpue)) %>%
  mutate(Source = ifelse(cpue_type == "cpue", "Observed", "Predicted")) %>% 
  mutate(month = lubridate::month(date),
         year = lubridate::year(date)) %>% 
  mutate(fad_free = month %in% c(7:9) & year > 2008) %>% 
  select(-cpue_type) %>% 
  # pivot_wider(names_from = Source, values_from = mean_cpue) %>% 
  group_by(year, Source) %>% 
  summarise(fad_effect = mean(mean_cpue[fad_free == TRUE]) / mean(mean_cpue[fad_free == FALSE]))


future_bet_fad_effect <- bet_model_comp_future$forest_predictions %>%
  mutate(cpue = (catch / days),
         date = lubridate::ymd(paste(year, month, "1", sep = '/')),
         mll = paste(month,lat,lon, sep = '-')) %>%
  filter(year > 2008,
         mll %in% consistent_bet_places$mll) %>%
  gather("cpue_type","cpue", cpue, cpue_hat) %>% 
  group_by(date, cpue_type) %>%
  summarise(mean_cpue = mean(cpue)) %>%
  mutate(Source = ifelse(cpue_type == "cpue", "Observed", "Predicted")) %>% 
  mutate(month = lubridate::month(date),
         year = lubridate::year(date)) %>% 
  mutate(fad_free = month %in% c(7:9) & year > 2008) %>% 
  select(-cpue_type) %>% 
  # pivot_wider(names_from = Source, values_from = mean_cpue) %>% 
  group_by(year, Source) %>% 
  summarise(fad_effect = mean(mean_cpue[fad_free == TRUE]) / mean(mean_cpue[fad_free == FALSE]))




true_bet_fad_effect_plot <- future_bet_fad_effect %>% 
  filter(Source == "Observed") %>% 
  ungroup() %>% 
  ggplot(aes(year, fad_effect)) + 
  geom_line()




bet_cpue_fit_plot <- bet_model_comp$forest_predictions %>%
  mutate(cpue = (catch / days),
         date = lubridate::ymd(paste(year, month, "1", sep = '/')),
         mll = paste(month,lat,lon, sep = '-')) %>%
  filter(year > 2008,
         mll %in% consistent_bet_places$mll) %>%
  gather("cpue_type","cpue", cpue, cpue_hat) %>% 
  group_by(date, cpue_type) %>%
  summarise(mean_cpue = mean(cpue)) %>%
  mutate(Source = ifelse(cpue_type == "cpue", "Observed", "Predicted")) %>%
  ungroup() %>%
  ggplot() +
  geom_rect(
    data = closed_frame,
    aes(
      xmin = date,
      xmax = date + 64,
      ymin = 0,
      ymax = 2.75
    ),
    fill = 'grey',
    alpha = 1
  ) +
  geom_line(aes(date, mean_cpue, color = Source)) +
  #facet_grid(method ~ .)+
  labs(x = "Date", y = "Mean CPUE (mt/day)", title = "A) Bigeye") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.title = element_blank())

future_bet_cpue_fit_plot <- bet_model_comp_future$forest_predictions %>%
  mutate(cpue = (catch / days),
         date = lubridate::ymd(paste(year, month, "1", sep = '/')),
         mll = paste(month,lat,lon, sep = '-')) %>%
  filter(year > 2008,
         mll %in% consistent_bet_places$mll) %>%
  gather("cpue_type","cpue", cpue, cpue_hat) %>% 
  group_by(date, cpue_type) %>%
  summarise(mean_cpue = mean(cpue)) %>%
  mutate(Source = ifelse(cpue_type == "cpue", "Observed", "Predicted")) %>%
  ungroup() %>%
  ggplot() +
  geom_rect(
    data = closed_frame,
    aes(
      xmin = date,
      xmax = date + 64,
      ymin = 0,
      ymax = 2.75
    ),
    fill = 'grey',
    alpha = 1
  ) +
  geom_line(aes(date, mean_cpue, color = Source)) +
  #facet_grid(method ~ .)+
  labs(x = "Date", y = "Mean CPUE (mt/day)", title = "A) Bigeye") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.title = element_blank())



# bet_cpue_fit_plot <- bet_model_comp$model_comp %>%
#   mutate(cpue = (catch / days),
#          date = lubridate::ymd(paste(year, month, "1", sep = '/'))) %>%
#   filter(method == "random-forest", year > 2008) %>%
#   group_by(date, catch_source, method) %>%
#   summarise(mean_cpue = mean(cpue)) %>%
#   mutate(Source = ifelse(catch_source == "catch", "Observed", "Predicted")) %>%
#   ungroup() %>%
#   ggplot() +
#   geom_rect(
#     data = closed_frame,
#     aes(
#       xmin = date,
#       xmax = date + 64,
#       ymin = 0,
#       ymax = 2.75
#     ),
#     fill = 'grey',
#     alpha = 1
#   ) +
#   geom_line(aes(date, mean_cpue, color = Source)) +
#   #facet_grid(method ~ .)+
#   labs(x = "Month", y = "Mean CPUE (mt/day)", title = "A) Bigeye") +
#   scale_y_continuous(expand = c(0, 0)) +
#   theme(legend.title = element_blank())

# SKJ

consistent_skj_places <- skj_dat %>% 
  filter(year < 2009,fad != "una") %>% 
  select(month,lat, lon) %>% 
  unique() %>% 
  mutate(mll = paste(month,lat,lon, sep = '-'))

skj_fad_effect <- skj_model_comp$forest_predictions %>%
  mutate(cpue = (catch / days),
         date = lubridate::ymd(paste(year, month, "1", sep = '/')),
         mll = paste(month,lat,lon, sep = '-')) %>%
  filter(year > 2008,
         mll %in% consistent_bet_places$mll) %>%
  gather("cpue_type","cpue", cpue, cpue_hat) %>% 
  group_by(date, cpue_type) %>%
  summarise(mean_cpue = mean(cpue)) %>%
  mutate(Source = ifelse(cpue_type == "cpue", "Observed", "Predicted")) %>% 
  mutate(month = lubridate::month(date),
         year = lubridate::year(date)) %>% 
  mutate(fad_free = month %in% c(7:9) & year > 2008) %>% 
  select(-cpue_type) %>% 
  # pivot_wider(names_from = Source, values_from = mean_cpue) %>% 
  group_by(year, Source) %>% 
  summarise(fad_effect = mean(mean_cpue[fad_free == TRUE]) / mean(mean_cpue[fad_free == FALSE]))


future_skj_fad_effect <- skj_model_comp_future$forest_predictions %>%
  mutate(cpue = (catch / days),
         date = lubridate::ymd(paste(year, month, "1", sep = '/')),
         mll = paste(month,lat,lon, sep = '-')) %>%
  filter(year > 2008,
         mll %in% consistent_bet_places$mll) %>%
  gather("cpue_type","cpue", cpue, cpue_hat) %>% 
  group_by(date, cpue_type) %>%
  summarise(mean_cpue = mean(cpue)) %>%
  mutate(Source = ifelse(cpue_type == "cpue", "Observed", "Predicted")) %>% 
  mutate(month = lubridate::month(date),
         year = lubridate::year(date)) %>% 
  mutate(fad_free = month %in% c(7:9) & year > 2008) %>% 
  select(-cpue_type) %>% 
  # pivot_wider(names_from = Source, values_from = mean_cpue) %>% 
  group_by(year, Source) %>% 
  summarise(fad_effect = mean(mean_cpue[fad_free == TRUE]) / mean(mean_cpue[fad_free == FALSE])) 

  true_skj_fad_effect_plot <- future_skj_fad_effect %>% 
  filter(Source == "Observed") %>% 
  ungroup() %>% 
  ggplot(aes(year, fad_effect)) + 
  geom_line()

future_skj_rmse <- future_skj_fad_effect %>% 
  ungroup() %>% 
  pivot_wider(names_from = Source, values_from = fad_effect) %>% 
  yardstick::mape(truth = Observed, estimate = Predicted) %>% 
  mutate(species = "Skipjack",
        holdout = "Post 2009 Data")


 skj_rmse <- skj_fad_effect %>% 
  ungroup() %>% 
  pivot_wider(names_from = Source, values_from = fad_effect) %>% 
  yardstick::mape(truth = Observed, estimate = Predicted) %>% 
   mutate(species = "Skipjack",
          holdout = "Summer FAD Closures")
 

 future_bet_rmse <- future_bet_fad_effect %>% 
   ungroup() %>% 
   pivot_wider(names_from = Source, values_from = fad_effect) %>% 
   yardstick::mape(truth = Observed, estimate = Predicted) %>% 
   mutate(species = "Bigeye",
          holdout = "Post 2009 Data")
 
 
 bet_rmse <- bet_fad_effect %>% 
   ungroup() %>% 
   pivot_wider(names_from = Source, values_from = fad_effect) %>% 
   yardstick::mape(truth = Observed, estimate = Predicted) %>% 
   mutate(species = "Bigeye",
          holdout = "Summer FAD Closures")
 
 mape_summary <- skj_rmse %>% 
   bind_rows(future_skj_rmse, future_bet_rmse, bet_rmse)
 
 fig_fad_effect_mape_plot <- mape_summary %>% 
   ggplot(aes(species, .estimate / 100, fill = holdout)) + 
   geom_col(position = "dodge") + 
   scale_fill_discrete(name = "Held-Out Data") + 
   scale_y_continuous(expand = expansion(c(0,0.05)), labels = percent, 
                      name = "Mean Absolute Percent Error") + 
   scale_x_discrete(name = '')
 
 
 ggsave(
   plot = fig_fad_effect_mape_plot,
   filename = file.path(run_dir, "fig_fad_effect_mape_plot.png"),
   height = 0.8 * fig_width,
   width = fig_width
 )
 
skj_cpue_fit_plot <- skj_model_comp$forest_predictions %>%
  mutate(cpue = (catch / days),
         date = lubridate::ymd(paste(year, month, "1", sep = '/')),
         mll = paste(month,lat,lon, sep = '-')) %>%
  filter(year > 2008,
         mll %in% consistent_skj_places$mll) %>%
  gather("cpue_type","cpue", cpue, cpue_hat) %>% 
  group_by(date, cpue_type) %>%
  summarise(mean_cpue = mean(cpue)) %>%
  mutate(Source = ifelse(cpue_type == "cpue", "Observed", "Predicted")) %>%
  ungroup() %>%
  ggplot() +
  geom_rect(
        data = closed_frame,
        aes(
          xmin = date,
          xmax = date + 64,
          ymin = 0,
          ymax = 30
        ),
        fill = 'grey',
        alpha = 1
      ) +
  geom_line(aes(date, mean_cpue, color = Source)) +
  #facet_grid(method ~ .)+
  labs(x = "Date", y = "Mean CPUE (mt/day)", title = "B) Skipjack") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.title = element_blank())


future_skj_cpue_fit_plot <- skj_model_comp_future$forest_predictions %>%
  mutate(cpue = (catch / days),
         date = lubridate::ymd(paste(year, month, "1", sep = '/')),
         mll = paste(month,lat,lon, sep = '-')) %>%
  filter(year > 2008,
         mll %in% consistent_skj_places$mll) %>%
  gather("cpue_type","cpue", cpue, cpue_hat) %>% 
  group_by(date, cpue_type) %>%
  summarise(mean_cpue = mean(cpue)) %>%
  mutate(Source = ifelse(cpue_type == "cpue", "Observed", "Predicted")) %>%
  ungroup() %>%
  ggplot() +
  geom_rect(
    data = closed_frame,
    aes(
      xmin = date,
      xmax = date + 64,
      ymin = 0,
      ymax = 30
    ),
    fill = 'grey',
    alpha = 1
  ) +
  geom_line(aes(date, mean_cpue, color = Source)) +
  #facet_grid(method ~ .)+
  labs(x = "Date", y = "Mean CPUE (mt/day)", title = "B) Skipjack") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.title = element_blank())


# Combine and export
fig_A1 <-
  ggpubr::ggarrange(
    bet_cpue_fit_plot + theme(axis.title.x = element_blank()),
    skj_cpue_fit_plot,
    ncol = 1,
    nrow = 2,
    common.legend = T,
    legend = "bottom",
    align = "v"
  )

ggsave(
  file = paste0(run_dir, "fig_A1", '.png'),
  fig_A1,
  height = 0.8 * fig_width,
  width = fig_width
)


fig_fa1 <-
  ggpubr::ggarrange(
    future_bet_cpue_fit_plot + theme(axis.title.x = element_blank()),
    future_skj_cpue_fit_plot,
    ncol = 1,
    nrow = 2,
    common.legend = T,
    legend = "bottom",
    align = "v"
  )

ggsave(
  file = paste0(run_dir, "fig_fa1", '.png'),
  fig_fa1,
  height = 0.8 * fig_width,
  width = fig_width
)



### Figure A2 - Predicted vs. Observed catches for both BET and SKJ (comparing models)
# BET
bet_fit_plot <- bet_model_comp$model_comp %>%
  group_by(year_month, catch_source, method) %>%
  summarise(catch = sum(catch)) %>%
  mutate(Source = ifelse(catch_source == "catch", "Observed", "Predicted")) %>%
  ggplot(aes(year_month, catch, color = Source)) +
  geom_line() +
  facet_grid(method ~ .) +
  labs(x = "Month", y = "Catch (mt)", title = "A) Bigeye") +
  theme(legend.title = element_blank()) +
  scale_y_continuous(labels = comma)

# SKJ
skj_fit_plot <- skj_model_comp$model_comp %>%
  group_by(year_month, catch_source, method) %>%
  summarise(catch = sum(catch)) %>%
  mutate(Source = ifelse(catch_source == "catch", "Observed", "Predicted")) %>%
  ggplot(aes(year_month, catch, color = Source)) +
  geom_line() +
  facet_grid(method ~ .) +
  labs(x = "Month", y = "Catch (mt)", title = "B) Skipjack") +
  theme(legend.title = element_blank()) +
  scale_y_continuous(labels = comma)

# Combine and export
fig_A2 <-
  ggpubr::ggarrange(
    bet_fit_plot + theme(axis.title.x = element_blank()),
    skj_fit_plot ,
    ncol = 1,
    nrow = 2,
    common.legend = T,
    legend = "bottom",
    align = "v"
  )

ggsave(
  file = paste0(run_dir, "fig_A2", '.png'),
  fig_A2,
  height = 1 * fig_width,
  width = fig_width
)

# Miscellaneous
# skj_cpue_fit_plot <- skj_model_comp$model_comp %>%
#   mutate(cpue = (catch / days),
#          date = lubridate::ymd(paste(year,month,"1",sep = '/'))) %>%
#   filter(method == "random-forest",year > 2008) %>%
#   group_by(date, catch_source, method) %>%
#   summarise(mean_cpue = mean((cpue))) %>%
#   spread(catch_source, mean_cpue) %>%
#   ungroup() %>%
#   # filter(method == 'random-forest') %>%
#   ggplot() +
#   geom_rect(data = closed_frame,
#             aes(
#               xmin = (date),
#               xmax = (date + 64),
#               ymin = 0,
#               ymax = 25),
#             fill = 'grey',
#             alpha = 0.05) +
#   geom_line(aes(date, catch_hat), size = 1) +
#   geom_point(aes(date, catch), size = 2, color = "red") +
#   labs(x = "Year", y = "Mean CPUE", title = "B) Skipjack")



# Model errors?
skj_model_comp$model_comp %>%
  group_by(method) %>%
  summarise(rmse = sqrt(mean(squared_error)))

bet_model_comp$model_comp %>%
  group_by(method) %>%
  summarise(rmse = sqrt(mean(squared_error)))


# Add in FAD Experiment ---------------------------------------------------

bet_reform <- bet_model_comp$forest_predictions

bet_fad_fishing <- bet_reform$fad != 'una' & bet_reform$fad != 'oth'

bet_no_fad <- bet_reform

bet_no_fad$fad[bet_fad_fishing] <-  'una'
#
# bet_reform <- bet_reform %>%
#   mutate(logit_pred_no_fad = predict(bet_model_comp$logit_forest, bet_no_fad, type = 'prob')[,'TRUE'],
#          pos_pred_no_fad =  predict(bet_model_comp$seen_forest, bet_no_fad),
#          log_cpue_hat_no_fad = logit_pred_no_fad * pos_pred_no_fad,
#          cpue_hat_no_fad = exp(log_cpue_hat_no_fad),
#          catch_hat_no_fad = cpue_hat_no_fad * days,
#          fad_effect = catch_hat_no_fad / catch_hat,
#          cpue_hat = exp(log_cpue_hat))

## NOTE GOT RID OF LOG TRANSFORM IN RF BUT DON'T WANT TO REWRITE NOW,
## HENCE ALL THE LOGS THAT AREN'T LOGS... SIGH.

bet_reform <- bet_reform %>%
  mutate(
    pos_pred_no_fad =  predict(
      bet_model_comp$seen_forest,
      bake(bet_model_comp$forest_recipe, bet_no_fad)
    )$.pred,
    log_cpue_hat_no_fad = pos_pred_no_fad,
    cpue_hat_no_fad = (log_cpue_hat_no_fad),
    catch_hat_no_fad = cpue_hat_no_fad * days,
    fad_effect = catch_hat_no_fad / catch_hat,
    cpue_hat = (log_cpue_hat)
  )


# bet_results <- bet_boot$fitted_model %>%
#   purrr::transpose() %>%
#   .$predictions %>%
#   bind_rows()


skj_reform <- skj_model_comp$forest_predictions

skj_fad_fishing <- skj_reform$fad != 'una' & skj_reform$fad != 'oth'

skj_no_fad <- skj_reform

skj_no_fad$fad[skj_fad_fishing] <-  'una'

# skj_reform <- skj_reform %>%
#   mutate(logit_pred_no_fad = predict(skj_model_comp$logit_forest, skj_no_fad, type = 'prob')[,'TRUE'],
#          pos_pred_no_fad =  predict(skj_model_comp$seen_forest, skj_no_fad),
#          log_cpue_hat_no_fad = logit_pred_no_fad * pos_pred_no_fad,
#          cpue_hat_no_fad = exp(log_cpue_hat_no_fad),
#          catch_hat_no_fad = cpue_hat_no_fad * days,
#          fad_effect = catch_hat_no_fad / catch_hat,
#          cpue_hat = exp(log_cpue_hat))


skj_reform <- skj_reform %>%
  mutate(
    pos_pred_no_fad =  predict(
      skj_model_comp$seen_forest,
      bake(skj_model_comp$forest_recipe, skj_no_fad)
    )$.pred,
    log_cpue_hat_no_fad = pos_pred_no_fad,
    cpue_hat_no_fad = (log_cpue_hat_no_fad),
    catch_hat_no_fad = cpue_hat_no_fad * days,
    fad_effect = catch_hat_no_fad / catch_hat,
    cpue_hat = (log_cpue_hat)
  )



# skj_reform <- skj_reform %>%
#   mutate(
#     fad_effect = cpue_hat_no_fad / cpue_hat,
#     no_fad_adjusted_true_cpue = cpue * fad_effect,
#     no_fad_catch = no_fad_adjusted_true_cpue * days
#   )
# 
# 
# bet_reform <- bet_reform %>%
#   mutate(
#     fad_effect = cpue_hat_no_fad / cpue_hat,
#     no_fad_adjusted_true_cpue = cpue * fad_effect,
#     no_fad_catch = no_fad_adjusted_true_cpue * days
#   )

skj_reform <- skj_reform %>%
  mutate(
    fad_effect = cpue_hat_no_fad / cpue_hat,
    no_fad_adjusted_true_cpue = cpue * fad_effect,
    no_fad_catch = cpue_hat_no_fad * days
  )


bet_reform <- bet_reform %>%
  mutate(
    fad_effect = cpue_hat_no_fad / cpue_hat,
    no_fad_adjusted_true_cpue = cpue * fad_effect,
    no_fad_catch = cpue_hat_no_fad * days
  )

bet_obs_v_pred_plot = bet_reform %>%
  ggplot(aes(catch, catch_hat)) +
  geom_point(shape = 21) +
  geom_smooth(method = 'lm') +
  geom_abline(aes(intercept = 0, slope = 1))

skj_obs_pred_plot = skj_reform %>%
  ggplot(aes(catch, catch_hat)) +
  geom_point(shape = 21) +
  geom_smooth(method = 'lm') +
  geom_abline(aes(intercept = 0, slope = 1))

# Flatten and collect predictions -----------------------------------------

flat_bet_reform =   bet_reform %>%
  select(year,
         month,
         lat,
         lon,
         fad,
         days,
         eez,
         pna_country,
         catch,
         no_fad_catch) %>%
  rename(bet_catch = catch, bet_no_fad_catch = no_fad_catch)

flat_skj_reform =   skj_reform %>%
  select(year,
         month,
         lat,
         lon,
         fad,
         days,
         eez,
         pna_country,
         catch,
         no_fad_catch) %>%
  rename(skj_catch = catch, skj_no_fad_catch = no_fad_catch)

flat_reform <- flat_bet_reform %>%
  left_join(
    flat_skj_reform,
    by = c(
      'year',
      'month',
      'lat',
      'lon',
      'fad',
      'days',
      'eez',
      'pna_country'
    )
  ) %>%
  mutate(
    skj_bought = skj_catch - skj_no_fad_catch,
    bet_saved = bet_catch - bet_no_fad_catch,
    bet_ratio = bet_saved - skj_bought
  ) %>% 
  filter(lat <= 25, lat >= -25)

bycatch_fiddle <- flat_reform %>%
  filter(year > 2000,
         fad != "una") %>%
  mutate(bycatch_rate = bet_catch / skj_catch)

bycatch_fiddle %>%
  summarise(
    med = median(bycatch_rate),
    min = min(bycatch_rate),
    max =  max(bycatch_rate)
  )

# Perform any filters


flat_reform <- flat_reform %>%
  filter(bet_saved > 0,
         skj_bought > 0,
         # year == max(year),!(month %in% drop_months),
         fad != 'oth') %>%  {
           if (pna_only == T) {
             filter(., pna_country == T & is.na(pna_country) == F)
             
           } else {
             .
           }
           
         }

model_predictions <- flat_reform

marginal_cost = flat_reform %>%
  mutate(cost = (skj_bought * skipjack_price) + bet_saved * bigeye_bycatch_price)

# Run Coaseian Bargain ----------------------------------------------------

# flat_reform <- flat_reform %>%
#   slice(1:10)

flat_reform <- flat_reform %>%
  dplyr::arrange(desc(bet_ratio))

fads_bought = 1:dim(flat_reform)[1]# seq(0, 1 * dim(flat_bet_reform)[1], length.out = 1000)


# generate biological fits ------------------------------------------------


if (run_biological_fits == TRUE) {
  species_names <- c("BET", "SKJ")
  
  fleet_names <- c("LL", "PL", "PS-FAD", "PS-UNA")
  
  results_dir <- run_dir
  
  # Apportion catches, effort, and selectivities for assessment year to the fleets we've defined above
  
  species_names <- species_names[1:2]
  
  fleet_dat <- FleetWrangling(
    species_names = species_names,
    data_year = assessment_year,
    fleet_names = fleet_names,
    results_dir = results_dir
  )
  
  fleet_dat$effort_fleet[, "LL"] <-
    fleet_dat$effort_fleet[, "LL"] * 0.001
  
  # Load assessment data into a life history object
  
  lh <- map(
    species_names,
    LoadLHParam,
    data_year = assessment_year,
    fleet_names = fleet_names,
    results_dir = results_dir
  )
  
  names(lh) <- species_names
  
  fitted_lh <- pmap(
    list(
      lh = lh,
      catch = fleet_dat$catch_fleet_species,
      selec = fleet_dat$selectivities
    ),
    TMBFit,
    data_year = assessment_year,
    fleet_names = fleet_names,
    effort = fleet_dat$effort_fleet,
    results_dir = results_dir,
    version = "msy",
    use_historic_recruits = FALSE
  )
  
  
  plot(fleet_dat$catch_fleet_species$BET[, 3],
       fleet_dat$effort_fleet[, 3])
  
  save(file = paste0(run_dir, "fit_data_", assessment_year, ".Rdata"),
       fitted_lh,
       fleet_dat)
  
  
  
} else {
  load(paste0(run_dir, "fit_data_", assessment_year, ".Rdata"))
  
  # load(here::here('data',paste0('fit_data_new_',assessment_year,'.Rdata')))
}


# process biological fits -------------------------------------------------


fitted_bet <- fitted_lh$BET

fitted_bet$ref_points %>%
  ggplot(aes(year, F_Fmsy, color = species)) +
  geom_point()


fitted_skj <- fitted_lh$SKJ

fitted_bet$fleet$selectivities$selectivity <-
  fleet_dat$selectivities$BET

run_time <- 200

max_age <- fitted_bet$lh$max_age

alpha <- fitted_bet$lh$alpha

beta <- fitted_bet$lh$beta

num_fleets <- 3

fitted_bet$catches_fleet <- fitted_bet$catches_fleet %>%
  as.data.frame() %>%
  mutate(year_quarter = rownames(.) %>% as.numeric()) %>%
  gather(fleet, catch,-year_quarter) %>%
  mutate(gear_type = fleet)


# fitted_bet$catches_fleet$gear_type <- fitted_bet$catches_fleet$fleet

fitted_bet$f_fleet_apical <- fitted_bet$f_fleet_apical %>%
  as.data.frame() %>%
  mutate(year_quarter = rownames(.) %>% as.numeric()) %>%
  gather(fleet, f,-year_quarter) %>%
  mutate(gear_type = fleet)

# fitted_bet$f_fleet_apical$gear_type <- fitted_bet$f_fleet_apical$fleet

# fitted_bet$fleet

fitted_bet$pop <- fitted_bet$pop %>%
  as.data.frame() %>%
  mutate(year_quarter = rownames(.) %>% as.numeric()) %>%
  gather(age, numbers,-year_quarter) %>%
  mutate(age = as.numeric(age))

fitted_bet$SSB_age <- fitted_bet$SSB_age %>%
  as.data.frame() %>%
  mutate(year_quarter = rownames(.) %>% as.numeric()) %>%
  gather(age, biomass,-year_quarter) %>%
  mutate(age = str_replace(age, '\\D', "") %>% as.numeric())



fleet_names <- unique(fitted_bet$catches_fleet$gear_type) %>%
  sort()
if (smooth_things == T) {
  fitted_bet$catches_fleet <- fitted_bet$catches_fleet %>%
    mutate(int_quarter = quarter - quarter %/% 1) %>%
    group_by(gear_type, int_quarter) %>%
    mutate(catch = RcppRoll::roll_mean(
      catch,
      align = 'right',
      n = 5,
      fill = NA,
      partial = T
    )) %>%
    arrange(gear_type) %>%
    ungroup() %>%
    select(-int_quarter)
  
  fitted_bet$f_fleet_apical <- fitted_bet$f_fleet_apical %>%
    mutate(int_quarter = quarter - quarter %/% 1) %>%
    group_by(gear_type, int_quarter) %>%
    mutate(f = RcppRoll::roll_mean(
      f,
      align = 'right',
      n = 5,
      fill = NA,
      partial = T
    )) %>%
    arrange(gear_type) %>%
    ungroup() %>%
    select(-int_quarter)
}

recent <- fitted_bet[c('catches_fleet',
                       'pop',
                       'SSB_age',
                       'f_fleet_apical')] %>%
  purrr::map(ungroup) %>%
  purrr::map( ~ mutate(., wtf = year_quarter)) %>%
  purrr::map( ~ mutate(., quarter = year_quarter)) %>%
  purrr::map( ~ mutate(., year = quarter %/% 1)) %>%
  purrr::map( ~ mutate(., int_quarter = (quarter - year + .25) / .25)) %>%
  purrr::map(~ filter(. , year == recent_year))


initial_n <-  recent$pop %>%
  filter(int_quarter == 4) %>%
  group_by(age) %>%
  summarise(numbers = mean(numbers)) %>%
  # select(-quarter, -year, - int_quarter, -year_quarter) %>%
  # select(-quarter, -year, - int_quarter, -year_quarter,-species) %>%
  spread(age, numbers) %>%
  as.matrix()

sel_at_age <- fleet_dat$selectivities$BET %>%
  as.data.frame() %>%
  mutate(fleet = rownames(.)) %>%
  gather(age, selectivity,-fleet) %>%
  mutate(age = as.integer(age)) %>%
  filter(fleet %in% unique(fitted_bet$catches_fleet$gear_type))

recent_recruits <- fitted_bet$recruitment_total %>%
  as.data.frame() %>%
  mutate(timestep = rownames(.) %>% as.numeric()) %>%
  mutate(year = floor(timestep)) %>%
  filter(year >= (max(year) - 5))

recent_recruits <- mean(recent_recruits$recruits)


fmsy <-
  tibble(fleet = names(fitted_bet$Fmsy_fleet),
         fmsy = fitted_bet$Fmsy_fleet)


fmsy <- recent$f_fleet_apical %>%
  group_by(fleet, gear_type, int_quarter) %>%
  summarise(f = mean(f)) %>%
  ungroup() %>%
  left_join(fmsy, by = "fleet") %>%
  mutate(f = fmsy) %>%
  select(-fmsy)

fmsy_run <- project_age_tuna(
  initial_n = initial_n,
  max_age = fitted_bet$max_age,
  alpha = fitted_bet$alpha,
  beta = fitted_bet$beta,
  f = fmsy,
  sel_at_age = sel_at_age,
  num_fleets = length(fleet_names),
  fleets = fleet_names,
  sp_at_age = fitted_bet$weight_at_age * fitted_bet$maturity_at_age,
  lh = fitted_bet,
  recruitment_form = "bh",
  recent_recruits = 0
)

fmsy_trajectory_plot <- fmsy_run$population %>%
  group_by(quarter) %>%
  summarise(ssb = sum(ssb), numbers = sum(numbers)) %>%
  mutate(recruits = fmsy_run$recruits %>% as.numeric()) %>%
  ggplot(aes(quarter, ssb / fitted_bet$SSBmsy)) +
  geom_line()


recent_f <- recent$f_fleet_apical %>%
  group_by(fleet, gear_type, int_quarter) %>%
  summarise(f = mean(f)) %>%
  ungroup()

status_quo <- project_age_tuna(
  initial_n = initial_n,
  max_age = fitted_bet$max_age,
  alpha = fitted_bet$alpha,
  beta = fitted_bet$beta,
  f = recent_f,
  sel_at_age = sel_at_age,
  num_fleets = length(fleet_names),
  fleets = fleet_names,
  sp_at_age = fitted_bet$weight_at_age * fitted_bet$maturity_at_age,
  lh = fitted_bet,
  run_time = run_time + 1,
  recruitment_form = "bh",
  recent_recruits = 0
)

ssb_trajectory_plot <- status_quo$population %>%
  group_by(quarter) %>%
  summarise(ssb = sum(ssb), numbers = sum(numbers)) %>%
  mutate(recruits = status_quo$recruits %>% as.numeric()) %>%
  ggplot(aes(quarter, ssb / fitted_bet$SSBmsy)) +
  geom_line()

ssb_trajectory_plot


sq_status <- status_quo$population %>%
  filter(quarter == max(quarter)) %>%
  summarise(ssb = sum(ssb)) %>%
  mutate(ssb_v_ssbmsy = ssb / fitted_bet$SSBmsy)


safe_coaseian_age_fad_bargain <- safely(coaseian_age_fad_bargain)
set.seed(42)
# mean_flat_reform <- flat_reform %>%
#   group_by(lat,lon,month) %>%
#   summarise()

mean_fad_effect <-  flat_reform %>%
  filter(year >= 2008) %>%
  group_by(fad, lat, lon, month) %>%
  summarise(
    skj_bought = mean(skj_bought),
    bet_saved = mean(bet_saved),
    skj_catch = mean(skj_catch),
    bet_catch = mean(bet_catch),
    days = mean(days),
    bet_no_fad_catch = mean(bet_no_fad_catch)
  ) %>%
  mutate(bet_ratio = bet_saved / skj_bought) %>%
  filter(bet_saved > 0,
         skj_bought > 0,!(month %in% drop_months),
         fad != 'oth',
         fad != 'una') %>%
  arrange(desc(bet_ratio)) %>%
  ungroup() %>%
  mutate(total_bet_saved = cumsum(bet_saved))

recentish_reported_catch <- recent$catches_fleet %>%
  filter(fleet == "PS-FAD") %>%
  summarise(tc = sum(catch)) %>% {
    .$tc
  }

mean_fad_effect <- mean_fad_effect %>%
  filter(total_bet_saved < recentish_reported_catch)


# run coaseian bargain ----------------------------------------------------



fads_bought <-  1:(dim(mean_fad_effect)[1])

if (fit_bargain == T) {
  if (run_in_parallel == F) {
    purchases <-  pblapply(
      fads_bought,
      coaseian_age_fad_bargain,
      dat = mean_fad_effect,
      initial_n = initial_n,
      alpha = fitted_bet$alpha,
      beta = fitted_bet$beta,
      num_fleets = length(fleet_names),
      fleet_names = fleet_names,
      f_by_fleet = recent_f,
      catch_by_fleet = recent$catches_fleet,
      sel_at_age = sel_at_age,
      numbers_at_age = recent$pop,
      lh = fitted_bet,
      bigeye_price = bigeye_price,
      bigeye_ps_price = bigeye_bycatch_price,
      sjk_price = skipjack_price,
      discount = disc_rate,
      time = run_time ,
      status_quo = status_quo,
      recruitment_form = "bh",
      include_bigeye = FALSE
    )
    
  } else{
    future::plan(future::multiprocess, workers = num_cores)
    
    purchases <- furrr::future_map(
      fads_bought,
      coaseian_age_fad_bargain,
      dat = mean_fad_effect,
      initial_n = initial_n,
      alpha = fitted_bet$alpha,
      beta = fitted_bet$beta,
      num_fleets = length(fleet_names),
      fleet_names = fleet_names,
      f_by_fleet = recent_f,
      catch_by_fleet = recent$catches_fleet,
      sel_at_age = sel_at_age,
      numbers_at_age = recent$pop,
      lh = fitted_bet,
      bigeye_price = bigeye_price,
      bigeye_ps_price = bigeye_bycatch_price,
      sjk_price = skipjack_price,
      discount = disc_rate,
      time = run_time,
      status_quo = status_quo,
      recruitment_form = "bh",
      include_bigeye = TRUE,
      .progress = TRUE
    )
  }
  
  purchases <- purchases %>%
    bind_rows()
  

  # save.image(file = paste(run_dir, 'coaseian_fad_bargain_results.Rdata', sep = ''))
  
  # save(file = paste(run_dir, 'coaseian_fad_bargain_results.Rdata', sep = ''),
  #      purchases)
  
  saveRDS(purchases,file = file.path(run_dir,"coaseian_fad_bargain_results.rds"))
  
} else{
  
  purchases <- readRDS(file = file.path(run_dir,"coaseian_fad_bargain_results.rds"))
  
}

purchases <-
  purrrlyr::dmap_at(purchases, which(str_detect(
    colnames(purchases), "npv|annuity|annualized"
  )), ~ .x / npv_factor, npv_factor = npv_factor)

purchases <- purchases %>%
  mutate(required_tax = pmax(0, (npv_costs - npv_benefits) / npv_benefits),
         tax_per_ton = bigeye_price * required_tax)

purchases$annuity_payments <- -pmin(0,purchases$npv_surplus) /  sum((1 + disc_rate)^(-(0:(run_time * 0.25 - 1))))
  

# Process Results ------------------------------------------------------------

purchases$id <- 1:nrow(purchases)

purchases <- purchases %>%
  filter(psf_f_message == "X-convergence (3)")

current_freeschool_skj_catch <- (
  skj_dat %>%
    group_by(year) %>%
    filter(fad == 'una') %>%
    summarise(catch  = sum(catch))  %>%
    filter(year == max(year)) %>%  {
      .$catch
    }
)

current_fad_skj_catch <- (
  skj_dat %>%
    group_by(year) %>%
    filter(fad != 'una') %>%
    summarise(catch  = sum(catch))  %>%
    filter(year == max(year)) %>%  {
      .$catch
    }
)

bau_skj_catch <- (skj_dat %>%
                    group_by(year) %>%
                    summarise(catch = sum(catch, na.rm = T)) %>%
                    filter(year == max(year)) %>%  {
                      .$catch
                    })

npv_bau_skipjack <-
  sum(bau_skj_catch * skipjack_price * (1 + disc_rate) ^ (-(0:run_time) * .25))

npv_bau_freeschool_skj <-
  sum(current_freeschool_skj_catch * skipjack_price * (1 + disc_rate) ^ (-(0:run_time) * .25))

npv_bau_fad_skj <-
  sum(current_fad_skj_catch * skipjack_price * (1 + disc_rate) ^ (-(0:run_time) * .25))


purchases = purchases %>%
  ungroup() %>%
  mutate(
    npv_new_freeschool = map_dbl(freeschool_skj,
                                 ~ sum(((current_freeschool_skj_catch + .x) * skipjack_price) * (1 + disc_rate) ^
                                         (-(0:run_time) * .25)
                                 )),
    annuity_per_fad_day = annuity_payments / fad_days_bght_per_year,
    canned_premium = ((npv_new_freeschool / npv_factor) - pmin(0, npv_surplus)) / (npv_new_freeschool / npv_factor) - 1
    # canned_premium = (-pmin(0, npv_surplus) - npv_new_fad) / npv_new_freeschool - 1
    # (npv_bau_skipjack - (npv_sq_skipjack + + (pmin(0, npv_surplus)))) /(npv_sq_skipjack + (pmin(0, npv_surplus)))
  )

# Process Data ------------------------ --------------------------------------

min_surp = signif(quantile(purchases$annualized_surplus, .1, na.rm = T), 1)

purchases <- purchases %>%
  as_tibble() %>%
  mutate(prop_fad_days = fad_days_bght_per_year / max(fad_days_bght_per_year)) %>%
  filter(npv_benefits > 0) %>%
  mutate(p_fad = fad_days_bght_per_year / max(fad_days_bght_per_year)) %>%
  mutate(
    ssb_v_ssb0 = ssb_final / fitted_bet$SSB0,
    ssb_v_ssbmsy = ssb_final / fitted_bet$SSBmsy
  )

# # # Height and width for combo maps (horizontal)
# combo_height = 6
# combo_width = 12
# 
# # # Height and width for facet figures (vertical)
# facet_height = 9.25
# facet_width = 9.25
# #
# # # Height and width for maps (square)
# map_height = 9.25
# map_width = 9.25
# #
# # # Height and width for special map
# smap_height = 6
# smap_width = 12

# Results Plots --------------------------------------------------------------

# -------------
### Figure 1
# WCPFC players
# -------------

# This plot is made in a different script because it's a mess and not dependant on this analysis


# --------------------------------
### Figure 2 
# Observed BET Bycatch rates (Map)
# --------------------------------

### Load shapefiles for map
pna_countries <- c("Palau",
                   "Papua New Guinea",
                   "Solomon Islands",
                   "Tuvalu",
                   "Nauru",
                   "Marshall Islands",
                   "Micronesia",
                   "Gilbert Islands",
                   "Line Group",
                   "Phoenix Group")

world <- map_data("world2")

map_layer <- geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey60", color = "grey40", size = 0.25)

eez <- sf::st_read(dsn = here::here("data", "world_eez_lr_-30_330"), layer = "world_eez_lr_-30_330")

eez_dat <- eez %>%  dplyr::filter(Territory1 %in% pna_countries)

eez_layer <- geom_sf(data = eez_dat, fill = NA, color = "grey40", size = 0.3)

## Observed PS Catch 
ps <-
  foreign::read.dbf(here("data", "catch_effort", "purse_seine_4", "PURSE_SEINE.DBF")) %>%
  purrr::set_names(tolower(colnames(.))) %>%
  as_tibble() %>%
  mutate(lat = as.numeric(str_extract(lat5, '[0-9]+')),
         lat_dir = str_extract(lat5, '[aA-zA]+'),
         lon = as.numeric(str_extract(lon5, '[0-9]+')),
         lon_dir = str_extract(lon5, '[aA-zA]+'))

# Convert to 0 - 360 degrees
ps$lat <- ifelse(ps$lat_dir == "S", ps$lat * -1, ps$lat)
ps$lon <- ifelse(ps$lon_dir == "W", (360 - ps$lon), ps$lon)


# # Extract effort and tidy
# ps_effort <- ps %>%
#   dplyr::select(yy, mm, lat, lon, days, contains("sets_")) %>%
#   gather(entry, effort, c(days, contains("sets_"))) %>%
#   mutate(effort_unit = str_replace(entry, "_.*", ""),
#          set_type = case_when(effort_unit == "sets" ~ str_replace(entry, ".*_", ""),
#                               TRUE ~ "NA")) %>%
#   dplyr::select(-entry)
# ps_effort$set_type[ps_effort$set_type == "NA"] <- NA

# Extract catch and tidy
ps_catch <- ps %>%
  dplyr::select(yy, mm, lat, lon, contains("_c_")) %>%
  gather(entry, catch, contains("_c_")) %>%
  mutate(species_set = str_replace(entry, "_c_", "_"),
         species = toupper(str_replace(species_set, "_.*", "")),
         set_type = str_replace(species_set, ".*_", "")) %>%
  dplyr::select(-entry, -species_set)

ps_catch_fad <- ps_catch %>%
  dplyr::filter(set_type %in% c("log", "dfad", "afad", "oth"))


# Get averages by lat/lon
ps_bycatch_reformat <- ps_catch_fad %>%
mutate(lat_center = lat + 2.5,
       lon_center = lon + 2.5) %>%
  dplyr::filter(yy >= 2000) %>%
  group_by(yy, mm, lat_center, lon_center, species) %>%
  summarize(catch = sum(catch, na.rm = T)) %>%
  dplyr::filter(species %in% c("BET", "SKJ")) %>%
  spread(species, catch)

ps_bycatch_spatial_median <- ps_bycatch_reformat %>%
  group_by(lat_center, lon_center) %>%
  summarize(bycatch_rate = median((BET) / (SKJ), na.rm = T)) %>%
  arrange(bycatch_rate) %>%
  na.omit()

ps_bycatch_monthly_median <- ps_bycatch_reformat %>%
  mutate(month = lubridate::month(mm, label = TRUE)) %>%
  group_by(month) %>%
  summarize(bycatch_rate = median((BET) / (SKJ), na.rm = T)) %>%
  arrange(bycatch_rate)

### Make spatial map
observed_bycatch_spatial_plot <- ggplot(ps_bycatch_spatial_median) +
  geom_tile(aes(x = lon_center, y = lat_center, fill = pmin(2, bycatch_rate)), alpha = 1) +
  eez_layer +
  scale_fill_viridis(option = "C",
                     guide = guide_colorbar(title = "Observed BET/SKJ catch", barwidth = 10, barheight = 1),
                     limits = c(0, max(max(ps_bycatch_spatial_median$bycatch_rate), max(ps_bycatch_monthly_median$bycatch_rate))))+
  map_layer+
  labs(x = "", y = "", title = "A)")+
  coord_sf(xlim = c(90, 220), ylim = c(-25, 25))




### Make monthly bar plot
observed_bycatch_monthly_plot <- ggplot(ps_bycatch_monthly_median) +
  geom_col(aes(x = month, y = bycatch_rate, fill = bycatch_rate))+
  scale_fill_viridis(option = "C",
                     guide = guide_colorbar(title = "Observed BET/SKJ catch", barwidth = 10, barheight = 1),
                     limits = c(0, max(max(ps_bycatch_spatial_median$bycatch_rate), max(ps_bycatch_monthly_median$bycatch_rate))))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = "", y = "Observed BET/SKJ catch", title = "B)")

fig_A3 <- ggarrange(observed_bycatch_spatial_plot, 
                   observed_bycatch_monthly_plot, 
                    ncol = 1, 
                    nrow = 2, 
                    align = "v", 
                    common.legend = T, 
                    legend = "bottom",
                    heights = c(5,2))

ggsave(
  file = paste0(run_dir, "fig_A3", '.png'),
  fig_A3,
  height = 1.1 * fig_width,
  width = fig_width
)

### Predicted PS bycatch rates 

# Coordinates represent SW corner of a 5x5 box - get centroids for plotting
bycatch_spatial_median <- flat_reform %>%
  mutate(long = ifelse(lon < 0, (360 + lon), lon),
         lat_center = lat + 2.5,
         lon_center = long + 2.5) %>%
  group_by(lat_center, lon_center) %>%
  summarize(bycatch_rate = median((bet_catch) / (skj_catch))) %>%
  arrange(bycatch_rate)

bycatch_monthly_median <- flat_reform %>%
  mutate(date = lubridate::ymd(paste(year, month, 1, sep = '/'))) %>%
  mutate(month = lubridate::month(date, label = TRUE)) %>%
  #filter(year > 2009) %>%
  group_by(month) %>%
  summarise( bycatch_rate = median(bet_catch / skj_catch))


### Make spatial map
bycatch_spatial_plot <- ggplot(bycatch_spatial_median) +
  geom_tile(aes(x = lon_center, y = lat_center, fill = pmin(2, bycatch_rate)), alpha = 1) +
  eez_layer +
  scale_fill_gradient(
    low = "white",
    high = "tomato",
    guide = guide_colorbar(
      barwidth = 15,
      barheight = 1,
      frame.colour = "black",
      title = expression(over("BET caught", "SKJ caught"))
    ),
    oob = squish,
    limits = c(0, max(
      max(bycatch_spatial_median$bycatch_rate),
      max(bycatch_spatial_median$bycatch_rate)
    ))
  ) + 
  map_layer + 
  theme(plot.margin = margin()) + 
  labs(x = "", y = "", title = "A)")+
  coord_sf(xlim = c(90, 220), ylim = c(-25, 25))
  

### Make monthly bar plot
bycatch_monthly_plot <- ggplot(bycatch_monthly_median) +
  geom_col(aes(x = month, y = bycatch_rate, fill = bycatch_rate))+
  scale_fill_gradient(
    low = "white",
    high = "tomato",
    guide = guide_colorbar(
      barwidth = 15,
      barheight = 1,
      frame.colour = "black",
      title = expression(over("BET caught", "SKJ caught"))
    ),
    oob = squish,
    limits = c(0, max(
      max(bycatch_spatial_median$bycatch_rate),
      max(bycatch_spatial_median$bycatch_rate)
    ))
  ) +
  scale_y_continuous(expand = c(0,0))+
  labs(x = "", y = "Median BET/SKJ catch", title = "B)") + 
  theme(plot.margin = margin())

fig_2 <- ggarrange(bycatch_spatial_plot, 
                    bycatch_monthly_plot, 
                    ncol = 1, 
                    nrow = 2, 
                    align = "v", 
                    common.legend = T, 
                    legend = "bottom",
                    heights = c(5,2))

ggsave(
  file = paste0(run_dir, "fig_2", '.png'),
  fig_2,
  height = 1.1 * fig_width,
  width = fig_width
)


# -------------------------------------------------------
### Figure 3
# Predicted BET saved/SKJ lost (Map and monthly bar plot)
# -------------------------------------------------------

purchases_spatial <- flat_reform %>%
  mutate(long = ifelse(lon < 0, (360 + lon), lon),
         lat_center = lat + 2.5,
         lon_center = long + 2.5) %>%
  mutate(bet_ratio_new = bet_saved / skj_bought,
         bycatch_rate = bet_catch / skj_catch) %>%
  arrange(bet_ratio_new)

purchases_spatial_median <- flat_reform %>%
  mutate(long = ifelse(lon < 0, (360 + lon), lon),
         lat_center = lat + 2.5,
         lon_center = long + 2.5) %>%
  group_by(lat_center, lon_center) %>%
  summarise(median_ratio_new = median((bet_saved) / (skj_bought)),
            bycatch_rate = median((bet_catch) / (skj_catch))) %>%
  arrange(median_ratio_new)

# Monthly value of BET saved/SKJ lost and bycatch rates
purchases_monthly_medians <- flat_reform %>%
  mutate(date = lubridate::ymd(paste(year, month, 1, sep = '/'))) %>%
  mutate(month = lubridate::month(date, label = TRUE)) %>%
  #filter(year > 2009) %>%
  group_by(month) %>%
  summarise(bet_ratio = median((bet_saved) / (skj_bought)),
            bycatch_rate = median(bet_catch / skj_catch))

### Make plots of log(BET saved/SKJ purchased)

# Lat/lon
purchase_map_plot <- ggplot(purchases_spatial_median) +
  geom_tile(aes(x = lon_center, y = lat_center, fill = median_ratio_new)) +
  eez_layer +
  coord_sf(xlim = c(90, 220), ylim = c(-25, 25)) +
  scale_fill_gradient2(
    guide = guide_colorbar(
      barwidth = 15,
      barheight = 1,
      frame.colour = "black",
      title = expression(over("BET Reduction", "SKJ Reduction"))
    ),
    oob = squish,
    limits = c(0, max(
      max(purchases_monthly_medians$bet_ratio),
      max(purchases_spatial_median$bet_ratio)
    ))
  ) +
  map_layer +
  labs(x = "", y = "", title = "A)")

# limits = c(0, max(max(purchases_monthly_medians$bet_ratio), max(purchases_spatial_median$bet_ratio)))


# Monthly
monthly_bet_value_plot <- purchases_monthly_medians %>%
  ggplot(aes(month, bet_ratio, fill = bet_ratio)) +
  geom_col() +
  scale_fill_gradient2(
                       guide = guide_colorbar(barwidth = 15,
                                              barheight = 1,
                                              title = expression(over("BET reduction", "SKJ reduction")),
                                              frame.colour = "black"),
                       limits = c(0, max(max(purchases_monthly_medians$bet_ratio), max(purchases_spatial_median$bet_ratio))),
                       oob = squish)+
  scale_y_continuous(expand = c(0,0))+
  labs(x = '', y = expression(over("BET saved", "SKJ lost")), title = "B)")


fig_3 <- ggarrange(purchase_map_plot, 
                   monthly_bet_value_plot, 
                    ncol = 1, 
                    nrow = 2, 
                    align = "v", 
                    common.legend = T, 
                    legend = "bottom",
                    heights = c(5,2))

ggsave(
  file = paste0(run_dir, "fig_3", '.png'),
  fig_3,
  height = 1.1 * fig_width,
  width = fig_width
)

### Figure 4 -------
#Four panel plot showing a) npv of surplus, b) BET tax, c) SKJ premium, d) conservation payments

fad_label = "% of FAD Days Purchased"

# A) NPV of surplus
# Bargain plot of surplus vs Fad days purchased
bargain_min <- quantile(purchases$npv_surplus)[4]
bargain_max <- round(max(purchases$npv_surplus) * 2)
range <-
  seq(
    from = plyr::round_any(bargain_min, 1, f = floor),
    to = plyr::round_any(bargain_max, 1, f = ceiling),
    by = 4
  )
range_labels <-
  c(paste("\u2264", dollar(range[1])), dollar(range[2:length(range)]))


max_ssb <- purchases$ssb_v_ssbmsy[which(purchases$npv_surplus <0)[1] - 1]

# Bargain plot of surplus vs SSB/SSBMSY as a function of FAD days purchased
conservation_plot <- purchases %>%
  mutate(ssb_v_ssbmsy = plyr::round_any(ssb_v_ssbmsy,.02)) %>% 
  group_by(ssb_v_ssbmsy) %>% 
  summarise(npv_surplus = mean(npv_surplus),
            p_fad = mean(p_fad)) %>% 
  ggplot(aes(ssb_v_ssbmsy, npv_surplus, fill = p_fad)) +
  geom_vline(aes(xintercept = max_ssb), color = "red", linetype = 2) + 
  geom_hline(aes(yintercept = 0), color = 'black', linetype = 1) +
  geom_col(width = 0.02) + 
  # geom_point(size = point_size) +
  # scale_y_continuous(breaks = range,
  #                    labels = range_labels,
  #                    limits = c(min(range), max(range)),
  #                    name = 'NPV of BET Surplus \n(US$, billions)', oob = squish) +
  scale_y_continuous(name = 'NPV of BET Surplus \n(US$, billions)', oob = squish, labels = dollar) +
  labs(x = expression(SSB / SSB[MSY])) +
  scale_fill_viridis(
    option = 'C',
    guide = guide_colorbar(
      title = fad_label,
      barwidth = 15,
      barheight = 1,
      frame.colour = "black"
    ),
    limits = c(0, 1),
    labels =  percent_format(accuracy = 1),
    direction = 1
  ) +
  labs(title = "A)")

# B) BET sashimi tax
tax_min <- min(purchases$required_tax)
tax_max <- max(purchases$required_tax)
tax_range <-
  seq(
    from = tax_min,
    to = plyr::round_any(tax_max, 1, f = ceiling),
    length.out = 4
  )

tax_plot <- purchases %>%
  mutate(ssb_v_ssbmsy = plyr::round_any(ssb_v_ssbmsy,.02)) %>% 
  group_by(ssb_v_ssbmsy) %>% 
  summarise(required_tax = mean(required_tax),
            p_fad = mean(p_fad)) %>% 
  ggplot(aes(ssb_v_ssbmsy, required_tax, fill = p_fad)) +
  geom_vline(aes(xintercept = max_ssb), color = "red", linetype = 2) + 
  geom_hline(aes(yintercept = 0), color = 'black', linetype = 1) +
  geom_col(width = 0.02) + 
  scale_x_continuous(name = expression(SSB / SSB[MSY])) +
  # scale_y_continuous(name = '% Tax on \nSashimi-Grade BET',
  #                    breaks = tax_range,
  #                    limits = c(min(tax_range), max(tax_range)),
  #                    labels = percent_format(tax_range, accuracy = 1)) +
  scale_y_continuous(name = '% Tax on \nSashimi-Grade BET',
                     labels = percent_format(accuracy = 1)) +
  scale_fill_viridis(
    option = 'C',
    guide = guide_colorbar(
      title = fad_label,
      barwidth = 15,
      barheight = 1,
      frame.colour = "black"
    ),
    limits = c(0, 1),
    labels =  percent_format(accuracy = 1),
    direction = 1
  ) +
  labs(title = "B)")

# C) SKJ canned premium
premium_min <- min(purchases$canned_premium)
premium_max <- max(purchases$canned_premium)
premium_range <-
  seq(
    from = premium_min,
    to = plyr::round_any(premium_max, 0.00000005, f = ceiling),
    length.out = 4
  )

canned_premium_plot <- purchases %>%
  mutate(ssb_v_ssbmsy = plyr::round_any(ssb_v_ssbmsy,.02)) %>% 
  group_by(ssb_v_ssbmsy) %>% 
  summarise(canned_premium = mean(canned_premium),
            p_fad = mean(p_fad)) %>% 
  ggplot(aes(ssb_v_ssbmsy, canned_premium, fill = p_fad)) +
  geom_vline(aes(xintercept = max_ssb), color = "red", linetype = 2) + 
  geom_hline(aes(yintercept = 0), color = 'black', linetype = 1) +
  geom_col(width = 0.02) + 
  scale_x_continuous(name = expression(SSB / SSB[MSY])) +
  scale_y_continuous(name = '% Premium on \nFree-School SKJ',
                     labels = percent_format(accuracy = 1)) +
  # scale_y_continuous(limits = c(min(premium_range), max(premium_range)),
  #                    breaks = premium_range,
  #                    labels = percent_format(premium_range, accuracy = 1),
  #                    name = '% Premium on \nFree-School SKJ') +
  scale_fill_viridis(
    option = 'C',
    guide = guide_colorbar(
      title = fad_label,
      barwidth = 15,
      barheight = 1,
      frame.colour = "black"
    ),
    limits = c(0, 1),
    labels =  percent_format(accuracy = 1),
    direction = 1
  ) +
  labs(title = "C)")

# D) Conservation payment plot vs SSB/SSBMSY as a function of FAD days purchased
conservation_payments_plot <- purchases %>%
  mutate(thing = pmax(0, annuity_per_fad_day * npv_factor)) %>% 
  mutate(ssb_v_ssbmsy = plyr::round_any(ssb_v_ssbmsy,.02)) %>% 
  group_by(ssb_v_ssbmsy) %>% 
  summarise(thing = mean(thing),
            p_fad = mean(p_fad)) %>% 
  ggplot(aes(ssb_v_ssbmsy, thing, fill = p_fad)) +
  geom_vline(aes(xintercept = max_ssb), color = "red", linetype = 2) + 
  geom_hline(aes(yintercept = 0), color = 'black', linetype = 1) +
  geom_col(width = 0.02) + 
  # geom_point(size = point_size) +
  scale_y_continuous(labels = dollar, name = 'Annuity Per FAD Day \n(US$)') +
  labs(x = expression(SSB / SSB[MSY])) +
  scale_fill_viridis(
    option = 'C',
    limits = c(0, 1),
    labels = percent_format(accuracy = 1),
    guide = guide_colorbar(
      title = fad_label,
      barwidth = 15,
      barheight = 1,
      frame.colour = "black"
    ),
    direction = 1
  ) +
  labs(title = "D)")

# Combine and export
fig_4 <-
  ggpubr::ggarrange(
    conservation_plot + theme(axis.title.x = element_blank()),
    tax_plot + theme(axis.title.x = element_blank()),
    canned_premium_plot + theme(axis.title.x = element_blank()),
    conservation_payments_plot,
    ncol = 2,
    nrow = 2,
    common.legend = T,
    legend = "bottom",
    align = "hv"
  )

ggsave(
  file = paste0(run_dir, "fig_4", '.png'),
  fig_4,
  height = 1.1 * fig_width,
  width = fig_width
)

### Figure 5 -----------
### NPV extinction plot
# Extinction Value Plot - If you assume that BAU leads to extinction, what's the value

ext_bargain_min <-
  floor(min(purchases$npv_extinction_surplus) / 100) * 100
ext_bargain_max <-
  ceiling(max(purchases$npv_extinction_surplus) / 100) * 100
ext_range <-
  seq(from = ext_bargain_min, to = ext_bargain_max, by = 2000)
# ext_range_labels <- c(paste("\u2264", dollar(range[1])), dollar(range[2:length(range)]))


fig_5 <- purchases %>%
  ggplot(aes(ssb_v_ssbmsy, npv_extinction_surplus, color = p_fad)) +
  geom_hline(aes(yintercept = 0), color = 'black',
             linetype = 1) +
  geom_vline(aes(xintercept = max(purchases$ssb_v_ssbmsy[purchases$npv_extinction_surplus >= 0])), color = "black", linetype = 2) +
  geom_point(size = point_size) +
  scale_x_continuous() +
  scale_y_continuous(name = 'NPV of Extinction Surplus \n(US$, billions)', oob = squish) +
  # scale_y_continuous(breaks = ext_range,
  #                    labels = dollar,
  #                    limits = c(ext_bargain_min, ext_bargain_max),
  #                    name = 'NPV of Extinction Surplus \n(US$, billions)', oob = squish) +
  labs(x = expression(SSB / SSB[MSY])) +
  scale_color_viridis(
    option = 'C',
    limits = c(0, 1),
    labels =  percent_format(accuracy = 1),
    guide = guide_colorbar(
      title = fad_label,
      barwidth = 15,
      barheight = 1,
      frame.colour = "black"
    ),
    direction = 1
  )

ggsave(
  file = paste0(run_dir, "fig_5", '.png'),
  fig_5,
  height = 0.6 * fig_width,
  width = fig_width
)

### Miscelaneous figures -------------

# Annuity plots
# Annual payment per FAD day vs fad days bought plot
new_annuity_vs_fads_plot <- purchases %>%
  ggplot(aes((fad_days_bght_per_year),
             (annuity_payments * 1e9) / fad_days_bght_per_year,
             color = p_fad
  )) +
  geom_point(size = point_size) +
  scale_y_continuous(name = 'Annual Payment per FAD Day (US$)', labels = dollar) +
  scale_x_continuous(name = 'FAD Days Purchased (days)',
                     labels = comma) +
  scale_color_viridis(
    option = 'C',
    guide = guide_colorbar(
      title = fad_label,
      barwidth = 15,
      barheight = 1,
      frame.colour = "black"
    ),
    limits = c(0, 1),
    labels =  percent_format(accuracy = 1),
    direction = 1
  )

# # Annual surplus vs payment per FAD day plot
# new_surplus_vs_annuity_plot <- purchases %>%
#   ggplot(aes(annuity_payments, annualized_surplus, color = p_fad)) +
#   annotate("rect", xmin = min(purchases$annuity_payments), xmax = max(purchases$annuity_payments[purchases$annualized_surplus >= 0]),
#            ymin = -Inf, ymax = Inf, alpha=0.3, fill="#95CE8AFF") +
#   geom_hline(aes(yintercept = 0), linetype = 1, color = "black") +
#   geom_vline(aes(xintercept = max(purchases$annuity_payments[purchases$annualized_surplus >= 0])), color = "#95CE8AFF", linetype = 2)+
#   geom_point(size = point_size, alpha = 0.5) +
#   scale_x_continuous(expand = c(0.01, 0), name = 'Annual Payments (million $US)',
#                      labels = dollar) +
#   scale_y_continuous(name = 'Annual Surplus (million $US)', labels = dollar) +
#   scale_color_viridis(option = 'C',
#                       guide = guide_colorbar(title = fad_label, barwidth = 15, barheight = 1),
#                       labels =  percent_format(accuracy = 1), direction = 1)

# Model plots
# bet_obs_v_pred_plot = bet_reform %>%
#   filter(log_catch > 1) %>%
#   ggplot(aes(catch, catch_hat)) +
#   geom_point(shape = 21) +
#   geom_smooth(method = 'lm', color = "#cc0052") +
#   labs(x = "Observed BET Catch (mt)", y = "Predicted BET Catch (mt)")+
#   scale_x_continuous(expand = c(0.01,0), limits = c(0, 2000))+
#   scale_y_continuous(expand = c(0.01,0), limits = c(0,2000))+
#   geom_abline(aes(intercept = 0, slope = 1))
#
# skj_obs_v_pred_plot = skj_reform %>%
#   # filter(log_catch > -2 & FAD == 1) %>%
#   ggplot(aes(catch, catch_hat)) +
#   geom_point(shape = 21) +
#   geom_smooth(method = 'lm', color = "#cc0052") +
#   labs(x = "Observed SKJ Catch (mt)", y = "Predicted SKJ Catch (mt)")+
#   scale_x_continuous(expand = c(0.01,0), limits = c(0, 30000))+
#   scale_y_continuous(expand = c(0.01,0), limits = c(0,30000))+
#   geom_abline(aes(intercept = 0, slope = 1))


### Management Actions vs Catch Plots -------------------------------------------------

## PS Catch and Effort Data
# ps <-
#   foreign::read.dbf(here("data", "catch_effort", "purse_seine_4", "PURSE_SEINE.DBF")) %>%
#   purrr::set_names(tolower(colnames(.))) %>%
#   as_tibble() %>%
#   rename(year = yy, month = mm)
# 
# ps_catch_ns <- ps %>%
#   tidyr::gather(catch_type, catch, contains('_c_')) %>%
#   select(-contains('sets')) %>%
#   mutate(
#     species = stringr::str_split(catch_type, '_c_', simplify = T)[, 1],
#     set_type = stringr::str_split(catch_type, '_c_', simplify = T)[, 2]
#   ) %>%
#   group_by(year, species) %>%
#   summarise(catch = sum(catch, na.rm = TRUE)) %>%
#   ungroup()
# 
# ps_effort_ns <- ps %>%
#   select(year, days, contains('sets')) %>%
#   group_by(year) %>%
#   summarize(effort = sum(sets_log + sets_dfad + sets_afad + sets_oth + sets_una)) %>%
#   arrange(year) %>%
#   mutate(gear = "PS")
# 
# ps_bet_ns <- ps_catch_ns %>%
#   filter(species == "bet") %>%
#   mutate("gear" = 'PS')
# 
# ## LL Catch and Effort Data
# 
# dbfer <- function(path) {
#   out <- foreign::read.dbf(path) %>%
#     purrr::set_names(colnames(.) %>% tolower()) %>%
#     dplyr::select(-dplyr::contains('_n')) %>%
#     gather(entry, catch, dplyr::contains('_c')) %>%
#     dplyr::mutate(species = toupper(str_replace(entry, '_c', '')),
#                   units = "mt") %>%
#     dplyr::select(-entry)
#   
#   return(out)
# }
# 
# paths <-
#   paste(
#     here("data", "catch_effort", "longline"),
#     "/LONGLINE_",
#     c('00', '60', '70', '80', '90'),
#     ".DBF",
#     sep = ''
#   )
# 
# # paths <- paste("./data/longline/LONGLINE_",c('00','60','70','80','90'),".DBF", sep = '')
# 
# ll <-  purrr::map_df(paths, dbfer) %>%
#   mutate(
#     lat = as.numeric(str_extract(lat5, '[0-9]+')),
#     lat_dir = str_extract(lat5, '[aA-zA]+'),
#     lon = as.numeric(str_extract(lon5, '[0-9]+')),
#     lon_dir = str_extract(lon5, '[aA-zA]+')
#   )
# 
# ll$lat <- ifelse(ll$lat_dir == "S", ll$lat * -1, ll$lat)
# ll$lon <- ifelse(ll$lon_dir == "W", (360 - ll$lon), ll$lon)
# 
# ll_effort <- ll %>%
#   rename(year = yy) %>%
#   dplyr::select(year, hhooks) %>%
#   group_by(year) %>%
#   summarize(effort = sum(hhooks) / 1000) %>%
#   arrange(year) %>%
#   mutate(gear = "LL")
# 
# ll_catch <- ll %>%
#   dplyr::select(yy, mm, lat, lon, catch, species, units) %>%
#   group_by(yy, species) %>%
#   summarise(catch = sum(catch, na.rm = TRUE)) %>%
#   ungroup()
# 
# bet_ll <- ll_catch %>%
#   filter(species == "BET") %>%
#   mutate("gear" = "LL")
# colnames(bet_ll) <- c("year", "species", "catch", "gear")
# 
# ## Combine LL and PS catches of BET and effort
# bet_catch <- rbind(ps_bet_ns, bet_ll) %>%
#   filter(year > 1995) %>%
#   select(-species) %>%
#   rename(fishery = gear, value = catch) %>%
#   arrange(desc(fishery)) %>%
#   mutate(type = "Catch")
# 
# effort <- rbind(ps_effort_ns, ll_effort) %>%
#   filter(year > 1995) %>%
#   rename(fishery = gear, value = effort) %>%
#   arrange(desc(fishery)) %>%
#   mutate(type = "Effort")
# 
# dat_combine <- bet_catch %>%
#   rbind(effort)
# 
# ## Load mgt data for both fisheries
# mgt <- read_csv(here("data", "mgt_data2.csv")) %>%
#   rename(fishery = Fishery,
#          mgt_name = mgt) %>%
#   mutate(xcen = st_yr + ((end_yr - st_yr) / 2),
#          ycen = ymin + ((ymax - ymin) / 2))
# 
# ## Do some cheating to make lines end nicely where I want them to
# bet_catch_subset <- bet_catch %>%
#   rename(lineend = value)
# 
# effort_subset <- effort  %>%
#   rename(lineend = value)
# 
# catch_fudge <-
#   data_frame(
#     year = c(2008.5, 2009.5),
#     type = "Catch",
#     fishery = "LL",
#     lineend = c(76679.49, 71191.84)
#   ) %>%
#   rbind(bet_catch_subset)
# 
# effort_fudge <-
#   data_frame(
#     year = c(2008.5),
#     type = "Effort",
#     fishery = "PS",
#     lineend = c(82138.93)
#   ) %>%
#   rbind(effort_subset)
# 
# mgt <- mgt %>%
#   left_join(catch_fudge,
#             by = c(
#               "st_yr" = "year",
#               "type" = "type",
#               "fishery" = "fishery"
#             )) %>%
#   left_join(effort_fudge,
#             by = c(
#               "st_yr" = "year",
#               "type" = "type",
#               "fishery" = "fishery"
#             ))
# mgt[is.na(mgt)] <- 0
# mgt <- mgt %>%
#   mutate(linee = lineend.x + lineend.y)
# 
# ## Split text so it will nicely fit in boxes based on \n in input
# mgt <- mgt %>%
#   mutate(
#     first = str_split(mgt_name, "\\\\n", simplify = T)[, 1],
#     second = str_split(mgt_name, "\\\\n", simplify = T)[, 2],
#     startbox = ifelse(fishery == "LL", 90000, 100000),
#     boxsize = ifelse(fishery == "LL", 10000, 20000),
#     textsize = 2000,
#     firsttext = ifelse(
#       second == "",
#       startbox + (boxsize * ycen),
#       ifelse(
#         fishery == "LL",
#         startbox + (boxsize * ycen) + 2000,
#         startbox + (boxsize * ycen) + 4000
#       )
#     ),
#     sectext = ifelse(
#       fishery == "LL",
#       startbox + (boxsize * ycen) - 2000,
#       startbox + (boxsize * ycen) - 4000
#     )
#   )
# 
# ## Actual plot
# a2 <- letters[seq(from = 1, to = length(unique(mgt$fishery)))]
# values2 <- c("Longline", "Purse Seine")
# facets2 <-
#   map2_chr(a2, values2, function(x, y) {
#     paste(x, ". ", y, sep = "")
#   })
# names(facets2) <- c("LL", "PS")
# 
# mgt_catch_facet_plot <- ggplot() +
#   geom_rect(
#     data = mgt,
#     aes(
#       xmin = st_yr,
#       xmax = end_yr,
#       ymin = startbox + (ymin * boxsize),
#       ymax = startbox + (ymax * boxsize),
#       fill = type
#     ),
#     color = "white",
#     alpha = 0.4
#   ) +
#   ggrepel::geom_text_repel(data = mgt,
#                            aes(x = xcen, y = firsttext, label = first),
#                            size = 3) +
#   ggrepel::geom_text_repel(data = mgt,
#                            aes(x = xcen, y = sectext, label = second),
#                            size = 3) +
#   geom_segment(data = mgt,
#                aes(
#                  x = st_yr,
#                  xend = st_yr,
#                  y = linee,
#                  yend = startbox + (ymax * boxsize),
#                  color = type,
#                  linetype = type
#                )) +
#   geom_line(data = dat_combine,
#             aes(x = year, y = value, color = type),
#             size = 1.5) +
#   geom_point(data = dat_combine,
#              aes(x = year, y = value, color = type),
#              size = 2) +
#   # geom_line(data = effort, aes(x = year, y = effort), size = 1.5, color = "#482677FF")+
#   # geom_point(data = effort, aes(x = year, y = effort), size = 2, color = "#482677FF")+
#   scale_x_continuous(lim = c(1996, 2019),
#                      breaks = seq(1996, 2018, by = 2)) +
#   scale_linetype_manual(values = c("dashed", "dotted")) + # first is "catch policies" second is "effort policies"
#   scale_y_continuous(
#     labels = comma,
#     sec.axis = sec_axis( ~ ., name = "Effort (10,000 hooks or sets)", labels = comma)
#   ) +
#   scale_color_manual(values = c("#1F968BFF", "#482677FF")) +
#   scale_fill_manual(values = c("#1F968BFF", "#482677FF")) +
#   labs(x = "Year", y = "Bigeye Catch (mt)") +
#   facet_wrap(
#     ~ fishery,
#     ncol = 1,
#     scales = "free_y",
#     labeller = as_labeller(facets2)
#   )



# make table --------------------------------------------------------------

purchases %>%
  mutate(rounded_s = round(ssb_v_ssbmsy,1)) %>%
  group_by(rounded_s) %>%
  summarise(Benefits = max(npv_benefits),
            Costs = max(npv_costs),
            Surplus = Benefits - Costs,
            `Percent of FADs Removed` = paste0(round(max(p_fad) * 100),'%')) %>%
  mutate(Benefits = paste0('$', prettyNum(Benefits, big.mark = ',')),
         Costs = paste0('$', prettyNum(Costs, big.mark = ',')),
         Surplus = paste0('$', prettyNum(Surplus, big.mark = ','))) %>%
  rename(`SSB/SSBmsy` = rounded_s) %>%
  knitr::kable(caption = 'NPV of revenue benefits (increases in bigeye catch) and revenue costs (losses in skipjack catch) required to achieve target levels of SSB/SSBmsy')

purchases %>%
  mutate(rounded_s = plyr::round_any(ssb_v_ssbmsy,.1, floor)) %>%
  group_by(rounded_s) %>%
  summarise(Benefits = signif(mean(npv_benefits),2),
            Costs = signif(mean(npv_costs),2),
            Surplus = Benefits - Costs,
            `Percent of FADs Removed` = paste0(round(max(p_fad) * 100),'%')) %>%
  mutate(Benefits = paste0('$', prettyNum(Benefits, big.mark = ',')),
         Costs = paste0('$', prettyNum(Costs, big.mark = ',')),
         Surplus = paste0('$', prettyNum(Surplus, big.mark = ','))) %>%
  rename(`SSB/SSBmsy` = rounded_s) %>%
  select(`Percent of FADs Removed`, everything()) %>% 
  write_csv(file.path(run_dir,"figure-table.csv")) 

### Data Wranging for Spatial Catch and Effort Maps ----------------------------------



# flat_reform %>%
#   mutate(date = lubridate::ymd(paste(year,month,1, sep = '/'))) %>%
#   mutate(month = lubridate::month(date, label = TRUE)) %>%
#   mutate(bet_ratio = pmin(2,(bet_saved) / (skj_bought))) %>%
#   ggplot(aes(x = bet_ratio, y=  factor(lat))) +
#   ggridges::geom_density_ridges()
#


### Save all miscellaneous plots ----------------------------------------------------------

plot_directory <- paste0(run_dir, "misc_plots/")
if (dir.exists(plot_directory) == F) {
  dir.create(plot_directory, recursive = T)
  
}

plot_list = ls()[grepl('_plot', ls())]

for (i in seq_along(plot_list)) {
  plot_plot <- get(plot_list[i])
  
  if (class(plot_plot)[1] == 'gg') {
    
      ggsave(
        file = paste(plot_directory, plot_list[i], '.png', sep = ''),
        plot_plot,
        height = fig_height,
        width = fig_width
      )

  } else{
    pdf(
      file = paste(plot_directory, plot_list[i], '.png', sep = ''),
      width = fig_width,
      height = fig_height
    )
    
    print(plot_plot)
    
    dev.off()
  }
  
}


# save plots

figures_to_save <- ls()[str_detect(ls(),"fig_")]

save(file = file.path(run_dir,"figures.Rdata"), list = figures_to_save)


rmarkdown::render(here::here("documents","figures-and-tables.Rmd"), params = list(run_name = run_name))

# save.image(file = paste(run_dir, 'coaseian_fad_bargain_figures.Rdata', sep = ''))
