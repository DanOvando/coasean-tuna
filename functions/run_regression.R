#' run_reg runs a specific regression on the tuna data
#'
#' @param dat raw data
#' @param dep_var dependent variable
#' @param logit_dep_var logit dependent variable
#' @param ind_vars independent variables
#' @param form type of regression ols tobit hurdle
#'
#' @return
#' @export

run_reg <-
  function(dat,
           dep_var,
           logit_dep_var = NA,
           ind_vars,
           form = 'ols') {
    reg_fmla <-
      as.formula(paste(dep_var, "~ ", paste(ind_vars, collapse = "+")))

    if (form == 'ols') {
      regression = RobustRegression(lm(reg_fmla, data = dat), dat = dat)
    }
    if (form == 'tobit') {
      regression = AER::tobit(reg_fmla, data = dat, left = min(dat[, dep_var]))
      # regression2 = vglm(reg_fmla, data = dat, tobit(Lower = min(dat[, dep_var])))

      regression$VCOV = vcov(regression)
    }
    if (form == 'glm'){
      regression <-  glm(
        reg_fmla ,
        family = gaussian(link = log),
        data = dat
      )

      regression$VCOV <-  vcov(regression)
    }

    return(regression)
  }
