
# species_names <- c("BET", "SKJ", "YFT")
# fleet_names <- c("LL", "PL", "PS-FAD", "PS-UNA")
# data_year <- 2014

### Function that wrangles catches and effort by species and fleet. 

FleetWrangling <- function(species_names,
                           data_year,
                           fleet_names,
                           results_dir){
  
  ### Setup --------------------------------
  num_sp <- length(species_names)
  num_fleets <- length(fleet_names)
  
  ## Paths to data
  dat_dir <- here::here("data", paste0("catch_effort/"))
  sp_dat_dir <- here::here("data", paste0(data_year, "_SA/"))
  
  ## Place to save output .csv files and data plots
  save_dir <- paste0(results_dir, data_year, "/", "inputs/")
  
  if (dir.exists(save_dir) == F) {
    dir.create(save_dir, recursive = T)
  }
  
  ## Quarterly cutoffs
  qtr_levels <- c(0,3,6,9,12) 
  qtr_labels <- c(1,2,3,4) 
  
  ### Create blank containers for catch/effort
  if(data_year == 2014){
    years <- c(1952, 2012)
  }else if(data_year == 2017){
    years <- c(1952, 2015)
  }
  
  year_quarter_seq <- seq(years[1], years[2]+0.75, by = 0.25)
  num_qtr <- length(year_quarter_seq)
  
  eff <- tibble(year_quarter = rep(year_quarter_seq, times = num_fleets),
                fleet = rep(fleet_names, each = num_qtr))
  
  cth <- tibble(year_quarter = rep(year_quarter_seq, times = num_fleets*num_sp),
                fleet = rep(fleet_names, each = num_qtr*num_sp),
                species = rep(species_names, each = num_qtr, length.out = num_fleets*num_sp*num_qtr))
  
  ### Assessment area boundaries ------------
  
  # Function to extract stock asessment regional boundaries by area
  AreaFun <- function(species){
    # Voodoo to get all 5x5 combinations of lat/lon in each region for each species
    out <- read_csv(paste0(sp_dat_dir, species, "/", species, "_regions_", data_year, ".csv")) %>% 
      group_by(id) %>%
      complete(nesting(id, region), lon = full_seq(lon, 5), lat = full_seq(lat, 5)) %>%
      ungroup() %>%
      group_by(region) %>%
      arrange(region, lon, lat) %>%
      dplyr::select(-id) %>%
      mutate(species = species)
  }
  
  areas <- map_df(species_names, AreaFun) 
  
  ### Gillnet (GN) fleet --- come back to this after we've dealt with ALB
  
  # if("GN" %in% fleet_names){
  # 
  # # All gillnet data
    # gn_dat <- foreign::read.dbf(paste0(dat_dir, "driftnet_0/DRIFTNET.DBF")) %>%
    #   set_names(colnames(.) %>% tolower()) %>%
    #   dplyr::select(-contains('_n')) %>%
    #   gather(entry, catch, contains('_c')) %>%
    #   mutate(species = toupper(str_replace(entry,'_c', ''))) %>%
    #   dplyr::select(-entry) %>%
    #   mutate(lat_num = as.numeric(str_replace(lat5,'[a-z, A-Z]', "")),
    #          lat_dir = str_extract(lat5, '[a-z, A-Z]+'),
    #          lon_num = as.numeric(str_replace(lon5, '[a-z, A-Z]', "")),
    #          lon_dir = str_extract(lon5, '[a-z, A-Z]+'),
    #          lat = ifelse(lat_dir == "N", lat_num, lat_num*-1),
    #          lon = ifelse(lon_dir == "E", lon_num, (180-lon_num)+180)) %>%
    #   dplyr::select(-lat_num, -lat_dir, -lon_num, -lon_dir) %>%
    #   mutate(quarter = as.numeric(cut(mm, qtr_levels, labels = qtr_labels)))
  #
  # # Gillnet effort
    # gn_effort <- gn_dat %>%
    #   group_by(yy, quarter, lat, lon, species) %>%
    #   summarize(effort = sum(days)) %>%
    #   ungroup() %>%
    #   group_by(yy, quarter, lat, lon) %>%
    #   summarize(effort = unique(effort))
  #
  # # Gillnet catches
  #
  # }
  
  ### Longline (LL) fleet -----------------------
  
  if("LL" %in% fleet_names){
    
    # Function to extract effort data from all longline .DBF files
    llEffortDbfer <- function(path){
      
      out <- foreign::read.dbf(path) %>% 
        set_names(colnames(.) %>% tolower()) %>%
        dplyr::select(-contains('_n'), -contains('_c')) %>%
        mutate(quarter = as.numeric(cut(mm, qtr_levels, labels = qtr_labels)),
               year_quarter = yy + (quarter/4 - 0.25))
      
      return(out)
    }
    
    # Function to extract catch data from all longline .DBF files
    llCatchDbfer <- function(path){
      
      out <- foreign::read.dbf(path) %>% 
        set_names(colnames(.) %>% tolower()) %>%
        dplyr::select(-contains('_n')) %>%
        gather(entry, catch, contains('_c')) %>% 
        mutate(species = toupper(str_replace(entry,'_c', ''))) %>%
        dplyr::select(-entry) %>%
        mutate(lat_num = as.numeric(str_replace(lat5,'[a-z, A-Z]', "")),
               lat_dir = str_extract(lat5, '[a-z, A-Z]+'),
               lon_num = as.numeric(str_replace(lon5, '[a-z, A-Z]', "")),
               lon_dir = str_extract(lon5, '[a-z, A-Z]+'),
               lat = ifelse(lat_dir == "N", lat_num, lat_num*-1),
               lon = ifelse(lon_dir == "E", lon_num, (180-lon_num)+180)) %>%
        dplyr::select(-lat_num, -lat_dir, -lon_num, -lon_dir) %>%
        mutate(quarter = as.numeric(cut(mm, qtr_levels, labels = qtr_labels)),
               year_quarter = yy + (quarter/4 - 0.25))
      
      return(out)
    }
    
    # Longline data path names
    ll_paths <- paste0(dat_dir, "longline/LONGLINE_", c('00', '60', '70', '80', '90'), ".DBF")
    
    # Longline effort
    ll_effort <- map_df(ll_paths, llEffortDbfer) %>%
      group_by(year_quarter) %>%
      summarize(effort = sum(hhooks)) %>% # LL EFFORT IS IN HOOKS
      mutate(fleet = "LL")
      
    # Longline catches
    ll_region <- map_df(ll_paths, llCatchDbfer) %>%
      group_by(year_quarter, lat, lon, species) %>%
      summarize(catch = sum(catch)) %>%
      dplyr::filter(species %in% species_names) %>%
      inner_join(areas, by = c("species", "lat", "lon"))
    
    ll_region_avg <- ll_region %>%
      group_by(species) %>%
      mutate(tot_catch = sum(catch)) %>%
      group_by(region, species) %>%
      summarize(region_catch = sum(catch),
                tot_catch = unique(tot_catch)) %>%
      ungroup() %>%
      mutate(region_prop = region_catch/tot_catch,
             fleet = "LL") %>%
      dplyr::select(region, fleet, species, region_prop)
    
    ll_catch <- ll_region %>%
      group_by(year_quarter, species) %>%
      summarize(catch = sum(catch)) %>%
      mutate(fleet = "LL")
    
  }
  
  ### Pole and Line (PL) fleet --------------------------
  
  if("PL" %in% fleet_names){
    
    # Pole and line effort
    pl_effort <- foreign::read.dbf(paste0(dat_dir, "pole_and_line_4/POLE_AND_LINE.DBF")) %>%
      set_names(colnames(.) %>% tolower()) %>%
      dplyr::select(-contains('_n'), -contains('_c')) %>%
      mutate(quarter = as.numeric(cut(mm, qtr_levels, labels = qtr_labels)),
             year_quarter = yy + (quarter/4 - 0.25)) %>%
      group_by(year_quarter) %>%
      summarize(effort = sum(days)) %>% # PL EFFORT IS IN DAYS
      mutate(fleet = "PL")
    
    # Pole and line catches
    pl_region <- foreign::read.dbf(paste0(dat_dir, "pole_and_line_4/POLE_AND_LINE.DBF")) %>%
      set_names(colnames(.) %>% tolower()) %>%
      dplyr::select(-contains('_n')) %>%
      gather(entry, catch, contains('_c')) %>%
      mutate(species = toupper(str_replace(entry,'_c', ''))) %>%
      dplyr::select(-entry) %>%
      mutate(lat_num = as.numeric(str_replace(lat5,'[a-z, A-Z]', "")),
             lat_dir = str_extract(lat5, '[a-z, A-Z]+'),
             lon_num = as.numeric(str_replace(lon5, '[a-z, A-Z]', "")),
             lon_dir = str_extract(lon5, '[a-z, A-Z]+'),
             lat = ifelse(lat_dir == "N", lat_num, lat_num*-1),
             lon = ifelse(lon_dir == "E", lon_num, (180-lon_num)+180)) %>%
      dplyr::select(-lat_num, -lat_dir, -lon_num, -lon_dir) %>%
      mutate(quarter = as.numeric(cut(mm, qtr_levels, labels = qtr_labels)),
             year_quarter = yy + (quarter/4 - 0.25)) %>%
      group_by(year_quarter, lat, lon, species) %>%
      summarize(catch = sum(catch)) %>%
      dplyr::filter(species %in% species_names) %>%
      inner_join(areas, by = c("species", "lat", "lon"))
    
    pl_region_avg <- pl_region %>%
      group_by(species) %>%
      mutate(tot_catch = sum(catch)) %>%
      group_by(region, species) %>%
      summarize(region_catch = sum(catch),
                tot_catch = unique(tot_catch)) %>%
      ungroup() %>%
      mutate(region_prop = region_catch/tot_catch,
             fleet = "PL") %>%
      dplyr::select(region, fleet, species, region_prop)
  
    pl_catch <- pl_region %>%
      group_by(year_quarter, species) %>%
      summarize(catch = sum(catch)) %>%
      mutate(fleet = "PL")
    
  }
  
  ### Purse seine (PS-UNA) fleet (unassociated) ----------------
  
  if("PS-UNA" %in% fleet_names){
    
    # Purse seine effort
    # We have to do some extra effort here because fishing days are not separated by set type - we therefore have to calculate the relative proportions of associated and unassociated sets and apply those proportions to days.
    
    ps_una_effort <- foreign::read.dbf(paste0(dat_dir, "purse_seine_4/PURSE_SEINE.DBF")) %>%
      set_names(colnames(.) %>% tolower()) %>%
      dplyr::select(-contains("_c_")) %>%
      group_by(yy, mm, lat5, lon5) %>%
      mutate(sets_fad = sets_log + sets_dfad + sets_afad + sets_oth,
             sets_tot = sets_una + sets_fad,
             prop_fad = ifelse(sets_tot == 0, 0, sets_fad/sets_tot),
             prop_una = ifelse(sets_tot == 0, 0, sets_una/sets_tot),
             days_una = round(prop_una * days, digits = 0),
             days_fad = round(prop_fad * days, digits = 0)) %>%
      ungroup() %>%
      dplyr::select(yy, mm, lat5, lon5, days_una, sets_una) %>%
      mutate(quarter = as.numeric(cut(mm, qtr_levels, labels = qtr_labels)),
             year_quarter = yy + (quarter/4 - 0.25)) %>%
      group_by(year_quarter) %>%
      summarize(effort = sum(days_una)) %>% # PS EFFORT IS IN DAYS
      mutate(fleet = "PS-UNA")
  
    # Purse seine catches
    ps_una_region <- foreign::read.dbf(paste0(dat_dir, "purse_seine_4/PURSE_SEINE.DBF")) %>%
      set_names(colnames(.) %>% tolower()) %>%
      dplyr::select(-contains("_log"), -contains("_dfad"), -contains("_afad"), -contains("_oth")) %>%
      gather(entry, catch, contains("_c_")) %>%
      mutate(species = toupper(str_replace(entry,'_c_.*', ''))) %>%
      dplyr::select(-entry) %>%
      mutate(lat_num = as.numeric(str_replace(lat5,'[a-z, A-Z]', "")),
             lat_dir = str_extract(lat5, '[a-z, A-Z]+'),
             lon_num = as.numeric(str_replace(lon5, '[a-z, A-Z]', "")),
             lon_dir = str_extract(lon5, '[a-z, A-Z]+'),
             lat = ifelse(lat_dir == "N", lat_num, lat_num*-1),
             lon = ifelse(lon_dir == "E", lon_num, (180-lon_num)+180)) %>%
      dplyr::select(-lat_num, -lat_dir, -lon_num, -lon_dir) %>%
      mutate(quarter = as.numeric(cut(mm, qtr_levels, labels = qtr_labels)),
             year_quarter = yy + (quarter/4 - 0.25)) %>%
      group_by(year_quarter, lat, lon, species) %>%
      summarize(catch = sum(catch)) %>%
      dplyr::filter(species %in% species_names) %>%
      inner_join(areas, by = c("species", "lat", "lon"))
    
    ps_una_region_avg <- ps_una_region %>%
      group_by(species) %>%
      mutate(tot_catch = sum(catch)) %>%
      group_by(region, species) %>%
      summarize(region_catch = sum(catch),
                tot_catch = unique(tot_catch)) %>%
      ungroup() %>%
      mutate(region_prop = region_catch/tot_catch,
             fleet = "PS-UNA") %>%
      dplyr::select(region, fleet, species, region_prop)
    
    ps_una_catch <- ps_una_region %>%
      group_by(year_quarter, species) %>%
      summarize(catch = sum(catch)) %>%
      mutate(fleet = "PS-UNA")
    
  }   
  
  ### Purse seine (PS-FAD) fleet (associated) ----------------
  
  if("PS-FAD" %in% fleet_names){
    
    # Purse seine effort
    # We have to do some extra effort here because fishing days are not separated by set type - we therefore have to calculate the relative proportions of associated and unassociated sets and apply those proportions to days.
    
    ps_fad_effort <- foreign::read.dbf(paste0(dat_dir, "purse_seine_4/PURSE_SEINE.DBF")) %>%
      set_names(colnames(.) %>% tolower()) %>%
      dplyr::select(-contains("_c_")) %>%
      group_by(yy, mm, lat5, lon5) %>%
      mutate(sets_fad = sets_log + sets_dfad + sets_afad + sets_oth,
             sets_tot = sets_una + sets_fad,
             prop_fad = ifelse(sets_tot == 0, 0, sets_fad/sets_tot),
             prop_una = ifelse(sets_tot == 0, 0, sets_una/sets_tot),
             days_una = round(prop_una * days, digits = 0),
             days_fad = round(prop_fad * days, digits = 0)) %>%
      ungroup() %>%
      dplyr::select(yy, mm, lat5, lon5, days_fad, sets_fad) %>%
      mutate(quarter = as.numeric(cut(mm, qtr_levels, labels = qtr_labels)),
             year_quarter = yy + (quarter/4 - 0.25)) %>%
      group_by(year_quarter) %>%
      summarize(effort = sum(days_fad)) %>% # PS EFFORT IS IN DAYS
      mutate(fleet = "PS-FAD")
    
    # Purse seine catches
    ps_fad_region <- foreign::read.dbf(paste0(dat_dir, "purse_seine_4/PURSE_SEINE.DBF")) %>%
      set_names(colnames(.) %>% tolower()) %>%
      dplyr::select(-contains("_una")) %>%
      gather(entry, catch, contains("_c_")) %>%
      mutate(species = toupper(str_replace(entry,'_c_.*', ''))) %>%
      dplyr::select(-entry) %>%
      mutate(lat_num = as.numeric(str_replace(lat5,'[a-z, A-Z]', "")),
             lat_dir = str_extract(lat5, '[a-z, A-Z]+'),
             lon_num = as.numeric(str_replace(lon5, '[a-z, A-Z]', "")),
             lon_dir = str_extract(lon5, '[a-z, A-Z]+'),
             lat = ifelse(lat_dir == "N", lat_num, lat_num*-1),
             lon = ifelse(lon_dir == "E", lon_num, (180-lon_num)+180)) %>%
      dplyr::select(-lat_num, -lat_dir, -lon_num, -lon_dir) %>%
      mutate(quarter = as.numeric(cut(mm, qtr_levels, labels = qtr_labels)),
             year_quarter = yy + (quarter/4 - 0.25)) %>%
      group_by(year_quarter, lat, lon, species) %>%
      summarize(catch = sum(catch)) %>%
      dplyr::filter(species %in% species_names) %>%
      inner_join(areas, by = c("species", "lat", "lon"))
    
    ps_fad_region_avg <- ps_fad_region %>%
      group_by(species) %>%
      mutate(tot_catch = sum(catch)) %>%
      group_by(region, species) %>%
      summarize(region_catch = sum(catch),
                tot_catch = unique(tot_catch)) %>%
      ungroup() %>%
      mutate(region_prop = region_catch/tot_catch,
             fleet = "PS-FAD") %>%
      dplyr::select(region, fleet, species, region_prop)
    
    ps_fad_catch <- ps_fad_region %>%
      group_by(year_quarter, species) %>%
      summarize(catch = sum(catch)) %>%
      mutate(fleet = "PS-FAD")
    
  }  
  
  ### Compile everything --------------------------------
  
  ## Effort
  effort_files <- ls()[grepl('_effort', ls())]

  effort_raw <- eval(as.name(effort_files[1]))
  
  for(i in 1:(length(effort_files)-1)) {
    
    effort_raw <- effort_raw %>%
      bind_rows(eval(as.name(effort_files[i+1])))
  }
  
  effort <- eff %>%
    left_join(effort_raw, by = c("year_quarter", "fleet")) %>%
    arrange(fleet, year_quarter)
  effort[is.na(effort)] <- 0
  
  # Save and plot effort by fleet
  write_csv(effort, paste0(save_dir, "effort_fleet.csv"))

  effort_plot_qtr <- ggplot(effort, aes(x = year_quarter, y = effort, color = fleet))+
    geom_point()+
    theme_bw()+
    theme(legend.position = "none")+
    scale_y_continuous(labels = comma)+
    facet_wrap(~fleet, ncol = 2, scales = "free_y")
  
  ggsave(paste0(save_dir, "effort_fleet_qtr.png"), effort_plot_qtr, dpi = 200)
  
  effort_plot_yr <- effort %>%
    mutate(year = floor(year_quarter)) %>%
    group_by(year, fleet) %>%
    summarize(effort = sum(effort)) %>%
    ggplot()+
    aes(x = year, y = effort, color = fleet)+
    geom_point()+
    theme_bw()+
    theme(legend.position = "none")+
    scale_y_continuous(labels = comma)+
    facet_wrap(~fleet, ncol = 2, scales = "free_y")
  
  ggsave(paste0(save_dir, "effort_fleet_yr.png"), effort_plot_yr, dpi = 200)
  
  # Convert to matrix
  effort_mat <- as.matrix(effort %>% spread(fleet, effort) %>% dplyr::select(-year_quarter))
  rownames(effort_mat) <- year_quarter_seq
  
  ## Catch 
  catch_files <- ls()[grepl('_catch', ls())]
  
  catch_raw <- eval(as.name(catch_files[1]))
  
  for(i in 1:(length(catch_files)-1)) {
    
    catch_raw <- catch_raw %>%
      bind_rows(eval(as.name(catch_files[i+1])))
  }
  
  catch <- cth %>%
    left_join(catch_raw, by = c("year_quarter", "species", "fleet")) %>%
    arrange(fleet, year_quarter, species)
  catch[is.na(catch)] <- 0
  
  # Export and make plots 
  write_csv(catch, paste0(save_dir, "catches_fleet_species.csv"))
  
  catch_plot_qtr <- ggplot(catch, aes(x = year_quarter, y = catch, fill = species))+
    geom_area(stat = "identity", position = "stack")+
    theme_bw()+
    scale_y_continuous(labels = comma)+
    facet_wrap(~fleet, ncol = 2, scales = "free_y")
  
  ggsave(paste0(save_dir, "catches_fleet_species_qtr.png"), catch_plot_qtr, dpi = 200)
  
  catch_plot_yr <- catch %>%
    mutate(year = floor(year_quarter)) %>%
    group_by(year, fleet, species) %>%
    summarize(catch = sum(catch)) %>%
    ggplot()+
    aes(x = year, y = catch, fill = species)+
    geom_area(stat = "identity", position = "stack")+
    theme_bw()+
    scale_y_continuous(labels = comma)+
    facet_wrap(~fleet, ncol = 2, scales = "free_y")
  
  ggsave(paste0(save_dir, "catches_fleet_species_yr.png"), catch_plot_yr, dpi = 200)
  
  # Function to extract catches by species and convert to matrices
  CatchMatrixExtract <- function(species_names, dat, qtrs){

    catch_mat <- as.matrix(dat %>% 
                             dplyr::filter(species == species_names) %>%
                             spread(fleet, catch) %>% 
                             dplyr::select(-year_quarter, -species))
    rownames(catch_mat) <- year_quarter_seq
    return(catch_mat)
    
  }
  
  catch_list <- map(species_names, CatchMatrixExtract, dat = catch, qtrs = year_quarter_seq)
  names(catch_list) <- species_names
  

  ### Selectivities -----------------
  
  # function to extract/combine selectivity files for all fleets
SelecExtract <- function(species, year, dat_dir){
  
  selec_dir <- paste0(dat_dir, species, "/", species, "_s_", year, "/")
  selec_files <- list.files(selec_dir)
  
  for(i in seq_along(selec_files)) {
    
    fleet <- str_extract(selec_files[i], ".+?(?=.txt)")
    
    if(!exists("selec")){
      selec <- read.table(paste0(selec_dir, selec_files[i]), header=TRUE, sep="\t") %>%
        mutate(SA_fleet = fleet,
               age_quarters = round(x))
    }
    
    else if (exists("selec")){
      temp_selec <- read.table(paste0(selec_dir, selec_files[i]), header=TRUE, sep="\t") %>%
        mutate(SA_fleet = fleet,
               age_quarters = round(x))
      selec <- rbind(selec, temp_selec)
      rm(temp_selec)
    }
  }
  
  selectivities <- selec %>%
    dplyr::select(-x) %>%
    rename(selectivity = y) %>%
    mutate(fleet = case_when(grepl("PS-UNA", SA_fleet) ~ "PS-UNA",
                             grepl("PS-ASS", SA_fleet) ~ "PS-FAD",
                             grepl("PS", SA_fleet) ~ "PS",
                             grepl("PL", SA_fleet) ~ "PL",
                             grepl("HL", SA_fleet) ~ "HL",
                             grepl("LL", SA_fleet) ~ "LL",
                             TRUE ~ "OTH"),
           region = as.numeric(str_extract(SA_fleet, "(\\d)+")),
           species = species)
  
  return(selectivities)
}

  selectivities <- map_df(species_names, SelecExtract, year = data_year, dat_dir = sp_dat_dir)
  
  selectivities_fleet_region <- selectivities %>%
    dplyr::filter(fleet %in% fleet_names) %>%
    group_by(region, fleet, species, age_quarters) %>%
    summarize(selectivity = mean(selectivity))
  
  # Weight selectivities by catches by region
  region_avg_files <- ls()[grepl('_region_avg', ls())]
  
  region_avg_raw <- eval(as.name(region_avg_files[1]))
  
  for(i in 1:(length(region_avg_files)-1)) {
    
    region_avg_raw <- region_avg_raw %>%
      bind_rows(eval(as.name(region_avg_files[i+1])))
  }
  
  region_avg <- region_avg_raw %>%
    arrange(region, fleet, species)
  
  selectivities_match <- selectivities_fleet_region %>%
    left_join(region_avg, by = c("region", "fleet", "species"))
  selectivities_match$region_prop[is.na(selectivities_match$region_prop)] <- 1
  
  # Calculate selectivities by fleet/species, weighted to the proportion of catches by that fleet/species coming from each region. 
  selectivities_fleet_weighted <- selectivities_match %>%
    group_by(fleet, species, age_quarters) %>%
    summarize(selec_raw = weighted.mean(selectivity, region_prop)) %>%
    mutate(selectivity = (selec_raw-min(selec_raw))/(max(selec_raw)-min(selec_raw))) # normalize so highest selectivity for each fleet/species combination is 1! 
    
  write_csv(selectivities_fleet_weighted, paste0(save_dir, "selectivities_fleet_species.csv"))
  
  selectivity_plot <- ggplot(selectivities_fleet_weighted, aes(x = age_quarters, y = selectivity, fill = species))+
    geom_bar(stat = "identity")+
    theme_bw()+
    facet_wrap(~ species + fleet, ncol = num_fleets, scales = "free_x")
  
  ggsave(paste0(save_dir, "selectivities_fleet_species.png"), selectivity_plot, dpi = 200)
  
  # Turn into a list of matrices by species
  SelectivityMatrixExtract <- function(species_names, dat, fleets){
    
    selec_mat <- as.matrix(dat %>% 
                             dplyr::select(-selec_raw) %>%
                             dplyr::filter(species == species_names) %>%
                             group_by(fleet) %>%
                             spread(age_quarters, selectivity) %>% 
                             ungroup() %>%
                             arrange(fleet) %>%
                             dplyr::select(-fleet, -species))
    rownames(selec_mat) <- fleets
    return(selec_mat)
    
  }
  
  selectivity_list <- map(species_names, SelectivityMatrixExtract, dat = selectivities_fleet_weighted, fleets = fleet_names)
  names(selectivity_list) <- species_names
  
 ### Return -----------------------------------
  
  return(list(effort_fleet = effort_mat,
              catch_fleet_species = catch_list,
              selectivities = selectivity_list))
  
}

