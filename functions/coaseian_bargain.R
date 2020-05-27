coaseian_bargain = function(skj_bought,
                            b_current,
                            f_current,
                            dat,
                            tunas,
                            c_current,
                            time = 50,
                            policy = 'constant f',
                            discount = 0.05,
                            bigeye_price = 10,
                            bigeye_ps_price = 1.5,
                            sjk_price = 1.5,
                            include_bigeye = F) {
  # skj_bought = 250
  # dat = bigeye_by_flag
  # tunas = tunas
  # time = 50
  # policy = 'constant f'

  dat$skj_bought = NA

  skj_to_buy = skj_bought
  for (i in 1:dim(dat)[1]) {
    dat$skj_bought[i] = min(dat$mean_targ_skipjack[i], skj_to_buy)

    skj_to_buy = max(0, skj_to_buy - dat$skj_bought[i])
  }

  dat = dat %>%
    mutate(bigeye_saved = skj_bought * mean_bycatch_rate)

  bigeye_saved = sum(dat$bigeye_saved, na.rm = T)

  status_quo = project_tuna(
    b_current = b_current,
    f_current = f_current,
    c_current = c_current,
    c_change = 0,
    g = tunas$g,
    phi = tunas$phi,
    time = time
  )

  projection = project_tuna(
    b_current = b_current,
    f_current = f_current,
    c_current = c_current,
    c_change = -bigeye_saved,
    g = tunas$g,
    phi = tunas$phi,
    time = time
  )

  bigeye_gains = pmax(0, projection$catch - status_quo$catch) # assume that all increases in catch are captured by the targeted fleet

  bigeye_benefits = bigeye_gains * bigeye_price

  if (include_bigeye == T) {
    skipjack_costs = c(0,
                       rep(skj_bought * sjk_price + bigeye_saved * bigeye_ps_price, time))
  } else{
    skipjack_costs = c(0,
                       rep(skj_bought * sjk_price, time))

  }

  npv_surplus = sum((bigeye_benefits - skipjack_costs) /  ((1 + discount) ^
                      (0:time)))

  out = data.frame(
    skj_bought = skj_bought,
    yield_surplus = last(bigeye_gains),
    npv_surplus = npv_surplus,
    b_final = last(projection$b_bmsy),
    f_final = last(projection$f_fmsy)
  )

  return(out)


}
