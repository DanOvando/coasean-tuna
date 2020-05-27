
####### General age structured population model ######

#' Function to grow age structured population. Fishing mortality can be specified, determined from catches, or determined from fleet dyanimcs model (in progress).
#'
#' \code{GrowStock} grows age structured population
#' @param purpose what the model is being used to calculate ('initial, 'Fmsy', 'EQ', predict catches', 'fit catches', or 'fleet model')
#' @param time time to run model (numeric value in number of quarters or 'EQ')
#' @param initial_pop population size at start of simulation (single value or vector of values equal to the number of age classes)
#' @param lh list of life history parameters
#' @param fleet list of fleet specific parameters (selectivities, catches, etc.)
#' @param fs fleet specific fishing mortalities (vector equal to the number of fleets)

GrowStock <-
  function(purpose,
           time,
           initial_pop,
           lh,
           fleet,
           fs_by_fleet = NA,
           tidy_results = F,
           catch_data = T) {
    ##### Set up run time #####
    
    if (time == 'EQ') {
      run_time <- 1000
      stop_time <- run_time + 1
    }
    
    if (time != 'EQ') {
      run_time <- time
      stop_time <- time + 1
    }
    
    ##### Set up storage for model outputs #####
    
    # Biological
    pop_storage <-
      matrix(0, nrow = stop_time, ncol = lh$max_age) #blank storage for population numbers. nrows = t+1, ncol = max age
    
    weight_storage <-
      matrix(0, nrow = stop_time, ncol = lh$max_age) #blank storage for population biomass. nrows = t+1, ncol = max age
    
    SSB_storage <-
      matrix(0, nrow = stop_time, ncol = 1) #blank storage for SSB.
    
    recruit_storage <-
      matrix(0, nrow = stop_time, ncol = 1) #blank storage for recruits
    
    # Fishing
    f_age_storage <-
      matrix(0, nrow = stop_time, ncol = lh$max_age) #blank storage for fishing mortality by age
    
    f_fleet_storage <-
      matrix(0, nrow = stop_time, ncol = length(fleet$names)) # blank storage for total fishing mortality.
    
    f_total_storage <-
      matrix(0, nrow = stop_time, ncol = 1) # blank storage for total fishing mortality.
    
    yield_age_storage <-
      matrix(0, nrow = stop_time, ncol = lh$max_age) #blank storage for fishing yields by age.
    
    yield_fleet_storage <-
      matrix(0, nrow = stop_time, ncol = length(fleet$names)) #blank storage for fishing yields by fleet.
    
    yield_total_storage <-
      matrix(0, nrow = stop_time, ncol = 1) #blank storage for total fishing yields.
    
    #
    # if(assessment == T){
    #
    #   # do something later
    #
    # }
    
    
    ##### Rename important parts of the lh files for easier use #####
    
    length_at_age <- lh$length_at_age
    weight_at_age <- lh$weight_at_age
    maturity_at_age <- lh$maturity_at_age
    
    ##### Merge fish and fleet data for later use #####
    
    merge_data <-
      MergeData(
        weight_storage[1, ],
        fleet,
        lh,
        t = 0,
        catch_data = catch_data,
        fs_by_fleet = c(0, 0, 0)
      )
    
    ##### Populate arrays with initial conditions #####
    
    if (purpose == "initial_EQ") {
      pop_storage[1, 1] <-
        lh$R0 # set first value of first row (population size of first age class for t = 0) to initial recruitment (first age class)
      
      recruit_storage[1] <- lh$R0
      
      pop_storage[1, ] <- lh$R0

        # pop_storage[1, 1] * exp(-lh$m * (0:(length(lh$m) - 1)))
      
       for (a in 2:(lh$max_age-1)){
      
         pop_storage[1,a] <- pmax(pop_storage[1,a-1]*(lh$s[a]), collapse_threshold) # population for each subsequent age class is the size of the previous times natural survivorship
       }
      
      pop_storage[1, lh$max_age] <-
        pmax(pop_storage[1, lh$max_age - 1] * ((lh$s[lh$max_age]) / (1 - (lh$s[lh$max_age]))), collapse_threshold) # plus group
      
      weight_storage[1, ] <-
        (pop_storage[1, ] * weight_at_age) / 1000 # converting number to biomass (and kg to mt)
      
      SSB_storage[1] <-
        sum(weight_storage[1, ] * maturity_at_age) # mature biomass
      
      f_fleet_storage[1, ] <- 0
      
      catch_results <-
        CalculateCatches(
          B = weight_storage[1, ],
          fleet = fleet,
          lh = lh,
          fs_by_fleet = f_fleet_storage[1, ],
          t = 0,
          merge_data = merge_data
        )
      f_age_storage[1, ] <- catch_results$f_by_age
      f_total_storage[1] <- sum(f_age_storage[1, ])
      
    } else{
      pop_storage[1, ] <-
        initial_pop #Seed population array with initial population
      
      recruit_storage[1] <- pop_storage[1, 1]
      
      weight_storage[1, ] <-
        (pop_storage[1, ] * weight_at_age) / 1000 # converting number to biomass (and kg to mt)
      
      SSB_storage[1] <-
        sum(weight_storage[1, ] * maturity_at_age) # mature biomass
      
      f_fleet_storage[1, ] <- fs_by_fleet
      
      catch_results <-
        CalculateCatches(
          B = weight_storage[1, ],
          fleet = fleet,
          lh = lh,
          fs_by_fleet = f_fleet_storage[1, ],
          t = 0,
          merge_data = merge_data
        )
      f_age_storage[1, ] <- catch_results$f_by_age
      f_total_storage[1] <- sum(f_age_storage[1, ])
      
    }
    
    
    ##### THIS IS THE MAIN GROWING/KILLING PART #####
    
    t <- 1 # start off time loop
    
    f_guess <- rep(lh$m, length = length(fleet_names))
    
    while (t < stop_time) {
      recruit_storage[t + 1] <-
        Recruitment_BH(SSB_storage[t], lh) #Calculate recruits w/ function
      
      merge_data <-
        MergeData(
          weight_storage[t, ],
          fleet,
          lh,
          t = t,
          catch_data = catch_data,
          fs_by_fleet = c(0, 0, 0)
        )
      
      ### Input catches or calculate catches by from selectivity and effort ###
      if (purpose == "initial_EQ") {
        
        f_fleet_storage[t + 1, ] <- 0
        
        catch_results <-
          CalculateCatches(
            B = weight_storage[t, ],
            fleet = fleet,
            lh = lh,
            fs_by_fleet = f_fleet_storage[t + 1, ],
            t = t,
            merge_data = merge_data
          )
        f_age_storage[t + 1, ] <- catch_results$f_by_age
        f_total_storage[t + 1] <- sum(f_age_storage[t + 1, ])
        
      } else if (purpose == "observed") {
        print(paste("Time is: ", t, " quarters", sep = ""))
        
        obs_fs_by_fleet <- nlminb(
          f_guess,
          Find_Fs_Fleet,
          lower = rep(0, length = length(fleet_names)),
          upper = rep(10, length = length(fleet_names)),
          B = weight_storage[t, ],
          fleet = fleet,
          lh = lh,
          t = t,
          merge_data = merge_data
        )

        obs_fs_by_fleet <- obs_fs_by_fleet$par
        
        f_guess <-  obs_fs_by_fleet
        
        f_fleet_storage[t + 1, ] <- obs_fs_by_fleet
        
        catch_results <-
          CalculateCatches(
            B = weight_storage[t, ],
            fleet = fleet,
            lh = lh,
            fs_by_fleet = obs_fs_by_fleet,
            t = t,
            merge_data = merge_data
          )
        
        f_age_storage[t + 1, ] <- catch_results$f_by_age
        f_total_storage[t + 1] <- sum(f_age_storage[t + 1, ])
        
      } else{
        print("Not a valid purpose")
      }
      
      ### Rename catch outputs for easier use ###
      
      survival_by_age <- catch_results$survival_by_age
      yield_by_fleet <- catch_results$yield_by_fleet
      yield_by_age <- catch_results$yield_by_age
      
      ### Take away catches to figure out how many fish we have left ###
      
      yield_age_storage[t + 1, ] <- yield_by_age
      
      yield_fleet_storage[t + 1, ] <- yield_by_fleet
      
      yield_total_storage[t + 1] <- sum(yield_age_storage[t + 1, ])
      
      pop_storage[t + 1, 1] <-
        recruit_storage[t + 1] #Add recruits to population
      
      plus_group <-
        pop_storage[t, lh$max_age] * survival_by_age[lh$max_age] #survival of fish already in plus group
      
      pop_storage[t + 1, 2:lh$max_age] <-
        pop_storage[t, 1:(lh$max_age - 1)] * survival_by_age[1:(lh$max_age - 1)] #Grow to the next age class
      
      pop_storage[t + 1, lh$max_age] <-
        pop_storage[t + 1, lh$max_age] + plus_group #add in fish that lived to enter plus group
      
      # pop_storage[t+1, lh$max_age] <- pop_storage[t, lh$max_age] * survival_by_age[lh$max_age] #Add in the survivors of the plus group
      
      weight_storage[t + 1, ] <-
        (pop_storage[t + 1, ] * weight_at_age) / 1000
      
      SSB_storage[t + 1] <-
        sum(weight_storage[t + 1, ] * maturity_at_age)
      
      ##### Stop things when population reaches equlibirum or the desired number of runs #####
      
      pop_change <- sum(pop_storage[t + 1, ]) - sum(pop_storage[t, ])
      
      #print(paste('Quarter', t, 'of', run_time))
      
      if (pop_change <= 0.00001 & t >= 200) {
        stop_time <- t + 1
        print(paste('Population reached EQ, time is:', t))
        
      } else if (pop_change > 0.00001 & t == 2000) {
        print(paste('Population failed to reach EQ, time is:', t))
        
      }
      
      t <- t + 1 #increase the time step by 1
      
    } # end time loop
    
    
    ##### Subset output files to reflect final stop time #####
    
    ### Biological results ###
    
    stock <- NULL
    
    stock$pop_by_age <- pop_storage[2:stop_time, ]
    stock$weight_by_age <- weight_storage[2:stop_time, ]
    stock$SSB <- SSB_storage[2:stop_time]
    stock$recruitment <- recruit_storage[2:stop_time]
    
    stock$final_num_at_age <- pop_storage[stop_time, ]
    stock$final_num <- sum(pop_storage[stop_time, ])
    stock$final_weight_at_age <- weight_storage[stop_time, ]
    stock$final_weight <- sum(weight_storage[stop_time, ])
    stock$final_SSB <- SSB_storage[stop_time]
    
    ### Performance results ###
    
    performance <- NULL
    
    performance$f_by_age <- f_age_storage[2:stop_time, ]
    performance$f_by_age_final <- f_age_storage[stop_time, ]
    performance$f_by_fleet <- f_fleet_storage[2:stop_time, ]
    performance$f_by_fleet_final <- f_fleet_storage[stop_time]
    performance$f_total <- f_total_storage[2:stop_time]
    performance$f_total_final <- f_total_storage[stop_time]
    
    performance$yields_by_age <- yield_age_storage[2:stop_time, ]
    performance$yields_by_fleet <- yield_fleet_storage[2:stop_time, ]
    performance$yields_total <- yield_total_storage[2:stop_time]
    
    if (tidy_results == T) {
      tidy_results <-
        TidyResults(
          stock = stock,
          performance = performance,
          lh = lh,
          fleet = fleet
        )
      
      return(list(
        stock = tidy_results$stock,
        performance = tidy_results$performance,
        lh = lh
      ))
      
    } else{
      return(list(
        stock = stock,
        performance = performance,
        lh = lh
      ))
    }
    
  } #Close population growth function
