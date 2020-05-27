find_bargain_f <-
  function(f_psf,
           f_by_fleet,
           numbers_at_age,
           lh,
           new_psfad_catch,
           fleet
  ){

    new_f <- f_by_fleet

    # f_psf_frame <- data_frame(int_quarter = 1:4, gear_type = 'PS-FAD', f_psf = f_psf)

    new_f$f[new_f$gear_type == fleet] <- f_psf

    new_f <- new_f %>%
      mutate(effective_f = f * selectivity)

    total_f_at_age <- new_f %>%
      ungroup() %>%
      mutate(effective_f = f * selectivity) %>%
      group_by(age,int_quarter) %>%
      dplyr::summarise(f = sum(effective_f))
    catch <- numbers_at_age  %>%
      left_join(total_f_at_age, by = c('age','int_quarter')) %>%
      mutate(catch = (f/(f + m)) * b_at_age * (1 - exp(-(f + m))))

    catch_by_fleet <-  new_f %>%
      left_join(catch %>% select(age,catch, int_quarter), by = c('age', 'int_quarter')) %>%
      group_by(age, int_quarter) %>%
      mutate(total_f_at_a = pmax(1e-6,sum(effective_f))) %>%
      ungroup() %>%
      mutate(catch_by_fleet = (effective_f / (total_f_at_a)) * catch) %>%
      group_by(gear_type,int_quarter) %>%
      dplyr::summarise(catch = sum(catch_by_fleet))

    psfad_catch <- catch_by_fleet$catch[catch_by_fleet$gear_type == fleet]

    obs_psfad_catch <- new_psfad_catch$new_catch[new_psfad_catch$gear_type == fleet]

    ss <- sum((log(psfad_catch + 1e-6) - log(obs_psfad_catch + 1e-6))^2)
    return(ss)

  }