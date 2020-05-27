


#' Age structured model for projecting tuna
#'
#' @param initial_n initial numbers at age
#' @param max_age max age
#' @param alpha alpha for BH recruit
#' @param beta beta for BH recruit
#' @param f current f by fleet
#' @param sel_at_age selectivity at age
#' @param num_fleets number of fleets
#' @param run_time number of quarters to run over
#' @param fleets names of fleets
#'
#' @return list of population trajectory
#' @export
#'
project_age_tuna <- function(initial_n,
                             max_age,
                             alpha,
                             beta,
                             f,
                             sel_at_age,
                             sp_at_age,
                             num_fleets,
                             run_time = 100,
                             fleets,
                             lh,
                             recruitment_form = "bh",
                             recent_recruits = NA)
{
  n_at_age <- matrix(0, run_time, max_age)

  ssb_at_age <- matrix(0, run_time, max_age)

  catch_at_age <- matrix(0, run_time, max_age)

  catch_at_fleet <- matrix(0, run_time, num_fleets)

  recruits <- vector(mode = 'double', length = run_time)

  n_at_age[1, ] <- initial_n

  ssb_at_age[1, ] <- initial_n * sp_at_age

  if (recruitment_form == "bh"){
  recruits[1] <-
    sum(ssb_at_age[1, ]) / (alpha + beta * sum(ssb_at_age[1, ]))
  } else {
    recruits <- rep(recent_recruits, run_time)
  }

  for (t in 2:run_time) {
    # if (t %% 4 >= 0) {
      n_at_age[t, 1] <- recruits[t - 1]
    # }

     last_int_quarter <- ifelse((t-1) %% 4 > 0,(t-1) %% 4 ,4)

     f_at_age <- sel_at_age %>%
       rename(gear_type = fleet) %>%
       left_join(
         f %>% filter(int_quarter == last_int_quarter) %>% dplyr::select(f, gear_type),
         by = c('gear_type')
       ) %>%
       mutate(f_at_age = as.numeric(f * selectivity)) %>%
       select(gear_type, age, f_at_age) %>%
       spread(age, f_at_age) %>%
       select(-gear_type) %>%
       as.matrix()

    total_f_at_age <- colSums(f_at_age)

    total_mort_at_age <- colSums(f_at_age) + lh$m

    survival_at_age <- exp(-total_mort_at_age)

    plus_group <-
      n_at_age[t - 1, max_age] * survival_at_age[max_age]

    n_at_age[t, 2:max_age] <-
      n_at_age[t - 1, 1:(max_age - 1)] * survival_at_age[1:(max_age - 1)]

    n_at_age[t, max_age] <-    n_at_age[t, max_age] + plus_group

    ssb_at_age[t, ] <-
      n_at_age[t, ] * lh$weight_at_age * lh$maturity_at_age

    catch_at_age[t - 1, ] <-
      lh$weight_at_age * (total_f_at_age / total_mort_at_age) *  n_at_age[t - 1,] * (1 - survival_at_age)

    catch_at_fleet[t - 1, ] <-
      rowSums(catch_at_age[t - 1, ] * f_at_age / total_f_at_age) %>% t()

    # if ((t %% 4) >= 0) {
    if (recruitment_form == "bh") {
      recruits[t] <-
        sum(ssb_at_age[t,]) / (alpha + beta * sum(ssb_at_age[t,]))
    }
    # }

  } #close time loop

  # if (t %% 4 >= 0) {
    n_at_age[t, 1] <- recruits[t - 1]
  # }

  last_int_quarter <- ifelse((t-1) %% 4 > 0,(t-1) %% 4 ,4)

  f_at_age <- sel_at_age %>%
    rename(gear_type = fleet) %>%
    left_join(f %>% filter(int_quarter == last_int_quarter) %>% dplyr::select(f, gear_type), by = c('gear_type')) %>%
    mutate(f_at_age = as.numeric(f*selectivity)) %>%
    select(gear_type,age,f_at_age) %>%
    spread(age, f_at_age) %>%
    select(-gear_type) %>%
    as.matrix()

  total_f_at_age <- colSums(f_at_age)

  total_mort_at_age <- colSums(f_at_age) + lh$m

  survival_at_age <- exp(-total_mort_at_age)

  catch_at_age[t, ] <-
    lh$weight_at_age * (total_f_at_age / total_mort_at_age) *  n_at_age[t,] * (1 - survival_at_age)

  catch_at_fleet[t, ] <-
    rowSums(catch_at_age[t, ] * f_at_age / total_f_at_age) %>% t()

  catch_at_fleet <-  catch_at_fleet %>%
    as_data_frame() %>%
    set_names(fleets) %>%
    mutate(quarter = 1:run_time) %>%
    tidyr::gather(fleet, catch, -quarter)

  catch_at_age <-  catch_at_age %>%
    as_data_frame() %>%
    set_names(1:max_age) %>%
    mutate(quarter = 1:run_time) %>%
    tidyr::gather(age, catch, -quarter)

  pop <- n_at_age %>%
    as_data_frame() %>%
    set_names(1:max_age) %>%
    mutate(quarter = 1:run_time) %>%
    tidyr::gather(age, numbers, -quarter)

  ssb <- ssb_at_age %>%
    as_data_frame() %>%
    set_names(1:max_age) %>%
    mutate(quarter = 1:run_time) %>%
    tidyr::gather(age, ssb, -quarter)

  pop <- pop %>%
    left_join(ssb, by = c('quarter', 'age'))

  out <- list(
    population = pop,
    recruits = recruits,
    catch_at_age = catch_at_age,
    catch_at_fleet = catch_at_fleet
  )
  return(out)

} #close function
