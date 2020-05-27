find_surplus_limits = function(i,
                               fad_effect,
                               flat_reform,
                               fads_bought,
                               tunas,
                               total_bigeye,
                               fleet_current,
                               bigeye_price,
                               bigeye_bycatch_price,
                               skipjack_price,
                               run_in_parallel = F,
                               alpha,
                               beta,
                               num_fleets = 3,
                               fleet_names,
                               f_by_fleet,
                               catch_by_fleet,
                               sel_at_age,
                               numbers_at_age,
                               lh,
                               fleet,
                               recent,
                               sjk_price,
                               discount) {

  flat_reform$skj_bought = flat_reform$true_skj_catch * fad_effect$skj_effect[i]

  flat_reform$bet_saved = flat_reform$true_bet_catch * fad_effect$bet_effect[i]

  if (run_in_parallel == F) {
    purchases = pblapply(
      fads_bought,
      coaseian_age_fad_bargain,
      dat = flat_reform,
      alpha = alpha,
      beta = beta,
      num_fleets = 3,
      fleet_names = fleet_names,
      f_by_fleet = recent$f_by_fleet,
      catch_by_fleet = recent$catch_by_fleet,
      sel_at_age = sel_at_age,
      numbers_at_age = recent$numbers_at_age ,
      lh = bet_lh,
      fleet = bet_fleet,
      recent = recent,
      bigeye_price = bigeye_price,
      bigeye_ps_price = bigeye_bycatch_price,
      sjk_price = skipjack_price,
      discount = disc_rate
    ) %>%
      bind_rows()
  } else{
    purchases = mclapply(
      fads_bought,
      coaseian_age_fad_bargain,
      dat = flat_reform,
      alpha = alpha,
      beta = beta,
      num_fleets = 3,
      fleet_names = fleet_names,
      f_by_fleet = recent$f_by_fleet,
      catch_by_fleet = recent$catch_by_fleet,
      sel_at_age = sel_at_age,
      numbers_at_age = recent$numbers_at_age ,
      lh = bet_lh,
      fleet = bet_fleet,
      recent = recent,
      bigeye_price = bigeye_price,
      bigeye_ps_price = bigeye_bycatch_price,
      sjk_price = skipjack_price,
      discount = disc_rate
    ) %>%
      bind_rows()
  }

  out = data.frame(fad_effect[i,], max_surplus = max(purchases$npv_surplus))

  return(out)

}