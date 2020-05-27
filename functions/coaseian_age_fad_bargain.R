coaseian_age_fad_bargain = function(fads_bought,
                                    dat,
                                    tunas,
                                    alpha,
                                    beta,
                                    initial_n,
                                    num_fleets = 3,
                                    fleet_names,
                                    f_by_fleet,
                                    catch_by_fleet,
                                    sel_at_age,
                                    numbers_at_age,
                                    lh,
                                    time = 50,
                                    discount = 0.05,
                                    bigeye_price = 10,
                                    bigeye_ps_price = 1.5,
                                    sjk_price = 1.5,
                                    include_bigeye = F,
                                    status_quo,
                                    recent_recruits = NA,
                                    recruitment_form = "historic") {
  # write(file = 'purchase_progress.txt', fads_bought, append = T)
  
  dat <- ungroup(dat)
  
  skj_bought = sum(dat$skj_bought[1:fads_bought], na.rm = T)
  
  freeschool_skj <-
    sum(dat$skj_catch[1:fads_bought], na.rm = T) - skj_bought #this is the worst name never. this is the new amount of FAF ASSOCIATED skipjack caught after moving that observation to freeschool
  
  # housekeeping step to enter how many bigeye will be "saved" from the FAD fishery
  bigeye_saved = dat %>%
    dplyr::slice(1:fads_bought) %>%
    mutate(int_quarter = ceiling(month / 3)) %>%
    group_by(int_quarter) %>%
    dplyr::summarise(bet_saved = sum(bet_saved)) %>%
    mutate(gear_type = 'PS-FAD')
  
  # housekeeping step to enter how many bigeye will be added to the freeschool fishery
  bigeye_still_caught <- dat %>%
    dplyr::slice(1:fads_bought) %>%
    mutate(int_quarter = ceiling(month / 3)) %>%
    group_by(int_quarter) %>%
    dplyr::summarise(bet_saved = -sum(bet_no_fad_catch)) %>%
    mutate(gear_type = 'PS-UNA')
  
  bigeye_change <- bigeye_saved %>%
    bind_rows(bigeye_still_caught)
  
  
  fad_days_bought = sum(dat$days[1:fads_bought], na.rm = T)
  
  # estimate the new catch be fleet as the difference in the currenta catch and the bet "saved". So, catch of fad fleet goes down, catch of free fleet goes up
  new_catch_by_fleet <- catch_by_fleet %>%
    left_join(bigeye_change, by = c('gear_type', 'int_quarter')) %>%
    mutate(new_catch =  pmax(0, catch - ifelse(is.na(bet_saved), 0 , bet_saved)))
  # Reduce F
  
  new_f <- f_by_fleet
  
  f_at_age_by_fleet <- sel_at_age %>%
    rename(gear_type = fleet) %>%
    left_join(new_f, by = 'gear_type') %>%
    filter(is.na(f) == F) %>%
    mutate(effective_f = f * selectivity)
  
  total_f_at_age <- f_at_age_by_fleet %>%
    group_by(age, int_quarter) %>%
    dplyr::summarise(f = sum(effective_f))
  
  quarter_n_at_age <- numbers_at_age  %>%
    left_join(data_frame(
      age = 1:length(lh$weight_at_age),
      weight_at_age = lh$weight_at_age,
      m = lh$m
    ),
    by = 'age') %>%
    mutate(b_at_age = numbers * weight_at_age)
  
  new_psf_f <- nlminb(
    new_f$f[new_f$gear_type == 'PS-FAD'] %>% jitter(),
    find_bargain_f,
    lower = 0,
    upper = 5,
    new_psfad_catch = new_catch_by_fleet,
    f_by_fleet = f_at_age_by_fleet,
    numbers_at_age = quarter_n_at_age,
    lh = lh,
    fleet = 'PS-FAD'
  )
  
  cc <- 0
  
  while (new_psf_f$convergence > 0 & cc < 9) {
    new_psf_f <- nlminb(
      new_f$f[new_f$gear_type == 'PS-FAD'] %>% jitter(),
      find_bargain_f,
      lower = 0,
      upper = 8,
      new_psfad_catch = new_catch_by_fleet,
      f_by_fleet = f_at_age_by_fleet,
      numbers_at_age = quarter_n_at_age,
      lh = lh,
      fleet = 'PS-FAD'
    )
    cc <-  cc + 1
    # print(cc)
    # if (cc > 10) {
    #
    #   browser()
    #
    #   stop("can't find f to justify new catch")
    # }
    
  }
  
  new_ps_f <- nlminb(
    new_f$f[new_f$gear_type == 'PS-UNA'],
    find_bargain_f,
    lower = 0,
    upper = 5,
    new_psfad_catch = new_catch_by_fleet,
    f_by_fleet = f_at_age_by_fleet,
    numbers_at_age = quarter_n_at_age,
    lh = lh,
    fleet = 'PS-UNA'
  )
  
  cc <-  0
  
  while (new_ps_f$convergence > 0 & cc < 9) {
    new_ps_f <- nlminb(
      new_f$f[new_f$gear_type == 'PS-UNA'],
      find_bargain_f,
      lower = 0,
      upper = 8,
      new_psfad_catch = new_catch_by_fleet,
      f_by_fleet = f_at_age_by_fleet,
      numbers_at_age = quarter_n_at_age,
      lh = lh,
      fleet = 'PS-UNA'
    )
    
    cc <- cc + 1
    # print(cc)
    
    # if (cc > 10) {
    #   browser()
    #   # stop("can't find f to justify new catch")
    # }
  }
  new_f$f[new_f$gear_type == 'PS-FAD'] <- new_psf_f$par
  
  new_f$f[new_f$gear_type == 'PS-UNA'] <- new_ps_f$par
  
  check_f <- sel_at_age %>%
    rename(gear_type = fleet) %>%
    left_join(new_f, by = 'gear_type') %>%
    filter(is.na(f) == F) %>%
    mutate(effective_f = f * selectivity) %>%
    ungroup()
  
  total_f <- sum(check_f$effective_f)
  #
  # hmm <- new_f
  #
  # new_f <- hmm
  #
  # new_f$f <-  new_f$f
  
  projection <- project_age_tuna(
    initial_n = initial_n,
    max_age = lh$max_age,
    alpha = alpha,
    beta = beta,
    f = new_f,
    sel_at_age = sel_at_age,
    num_fleets = num_fleets,
    fleets = fleet_names,
    sp_at_age = lh$weight_at_age * lh$maturity_at_age,
    lh = lh,
    run_time = time + 1,
    recent_recruits = recent_recruits,
    recruitment_form = recruitment_form
  )
  
  bigeye_gains <-
    
    (
      projection$catch_at_fleet %>% filter(fleet == 'LL') %>% select(catch) -
        status_quo$catch_at_fleet %>% filter(fleet == 'LL') %>% select(catch)
    )
  
  bigeye_gains <- pmax(0, bigeye_gains$catch) %>% as.numeric()
  
  # assume that all increases in catch are captured by the targeted fleet
  
  bigeye_benefits <-  bigeye_gains * bigeye_price
  
  extinction_revenue <-  (projection$catch_at_fleet %>%
                            filter(fleet == 'LL') %>%
                            select(catch) %>%
                            {
                              .$catch
                            }) * bigeye_price
  
  # benefits are quarterly but costs are annual, so divide costs by 4
  if (include_bigeye == T) {
    skipjack_costs =  rep(skj_bought * sjk_price + sum(bigeye_saved$bet_saved) * bigeye_ps_price,
                          time + 1) * 0.25
  } else{
    skipjack_costs = rep(skj_bought * sjk_price, time + 1) * 0.25
  }
  
  discount_factor <- (1 + discount) ^ (-seq(0, run_time) * 0.25)
  
  npv_benefits <-  sum(bigeye_benefits * discount_factor)
  
  npv_extinction_revenue <-
    sum(extinction_revenue * discount_factor)
  
  npv_costs <-  sum(skipjack_costs * discount_factor)
  
  npv_surplus <-
    sum((bigeye_benefits - skipjack_costs) * discount_factor)
  
  npv_extinction_surplus <- npv_extinction_revenue - npv_costs
  
  annualized_payments	<- -mean(pmin(0, npv_benefits - npv_costs))
  # annualized_payments	<- mean(skipjack_costs * (1 + discount) ^ -(time * 0.25))
  # (sum(skipjack_costs * discount_factor) * discount) / (1 - (1 + discount) ^
  #                                                         -(time * .25))
  
  annualized_surplus	<-
    (npv_surplus * discount) / (1 - (1 + discount) ^
                                  -(time * .25))
  
  final_ssb <-  projection$population %>%
    filter(quarter == max(quarter)) %>%
    dplyr::summarise(ssb = sum(ssb))
  
  final_catch <-
    projection$catch_at_age %>%
    filter(quarter == max(quarter)) %>%
    dplyr::summarise(catch = sum(catch))
  
  out <-  data_frame(
    skj_bought = skj_bought,
    freeschool_skj = freeschool_skj,
    bet_saved = sum(bigeye_saved$bet_saved),
    fad_days_bght_per_year = fad_days_bought,
    annuity_payments = annualized_payments,
    annualized_surplus = annualized_surplus,
    yield_surplus = last(bigeye_gains),
    npv_surplus = npv_surplus,
    npv_benefits = npv_benefits,
    npv_costs = npv_costs,
    npv_extinction_revenue = npv_extinction_revenue,
    npv_extinction_surplus = npv_extinction_surplus,
    ssb_final = final_ssb$ssb,
    catch_final = final_catch$catch,
    total_f = total_f,
    psf_f_objective_function = new_psf_f$objective,
    ps_f_objective_function = new_ps_f$objective,
    psf_f_message = new_psf_f$message,
    ps_f_message = new_ps_f$message
  )
  # if (fads_bought == 26){
  #
  #   browser()
  # }
  return(out)
  
}
