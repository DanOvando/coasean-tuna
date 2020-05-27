#' project_tuna runs a pella-tomlinson model forward using a given policy and starting conditions
#'
#' @param b_current current B/Bmsy
#' @param f_current current F/Fmsy
#' @param c_current current catch
#' @param c_change proposed change in catch
#' @param g pt growth parameter
#' @param phi pt production parameter
#' @param policy policy to be used
#' @param time number of years to run policy for
#' @param current_year the value of the current year
#'
#' @return a data.frame with projected tuna population and catches
#' @export
#'
project_tuna = function(b_current,
                        f_current,
                        c_current,
                        fleet_current = NA,
                        c_change = 9,
                        g = 0.1,
                        phi = 0.188,
                        policy = 'constant f',
                        time = 50,
                        current_year = 2014) {
  tuna_frame = data.frame(
    year = current_year:(current_year + time),
    b_bmsy = NA,
    f_fmsy = NA,
    c_msy = NA,
    catch = NA,
    targ_catch = NA
  )

  tuna_frame$b_bmsy[1] = b_current

  tuna_frame$f_fmsy[1] = f_current

  tuna_frame$catch = c_current

  msy = c_current / (b_current * f_current)

  tuna_frame$c_msy[1] = c_current / msy

  fleet_f = (fleet_current / msy) / b_current

  prop_targeted = fleet_f$bycatch / sum(fleet_f)

  tuna_frame$targ_catch[1] = c_current * prop_targeted

  if (policy == 'constant f') {
    new_f = ((c_current + c_change) / msy) / b_current

    new_fleet_f = ((fleet_current + c(c_change,0)) / msy) /b_current

    prop_targeted = new_fleet_f$bycatch / sum(new_fleet_f)

    tuna_frame$f_fmsy[2:(time + 1)] = new_f
  }

  for (t in 2:(time + 1)) {
    tuna_frame$b_bmsy[t] =  tuna_frame$b_bmsy[t - 1] + ((phi + 1) / phi) * (g *  tuna_frame$b_bmsy[t - 1]) * (1 - ((tuna_frame$b_bmsy[t - 1] ^ phi) / (phi + 1))) - g * tuna_frame$b_bmsy[t - 1] * tuna_frame$f_fmsy[t - 1]

    if (policy == 'constant catch') {
      new_f = ((c_current + c_change) / msy) / tuna_frame$b_bmsy[t]

      tuna_frame$f_fmsy[t] = new_f

    }
    tuna_frame$c_msy[t] = tuna_frame$b_bmsy[t] * tuna_frame$f_fmsy[t]

    tuna_frame$catch[t] = tuna_frame$c_msy[t] * msy

    tuna_frame$targ_catch[t] = tuna_frame$catch[t] * prop_targeted
  }

  tuna_frame$msy = msy

  return(tuna_frame)

}
