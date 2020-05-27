#' fit_fad runs a regression and predicts changes in skj and bet catch from fad reduction
#'
#' @param i index of row to be predicted
#' @param reg_dat regression data
#' @param dep_var dependent variable
#' @param ind_vars independent variables
#' @param form type of regression to use ols tobit hurdle
#'
#' @return
#' @export
#'
fit_fad <- function(i,
                    reg_dat,
                    dep_var,
                    ind_vars,
                    form = 'ols',
                    gam_smooth) {
  if (is.character(i)) {
    training = reg_dat

    test = reg_dat
  } else{
    training = reg_dat[-i,]

    test = reg_dat[i,]
  }

  if (form != 'hurdle') {
    reg = run_reg(
      dat = training,
      dep_var = dep_var,
      ind_vars = ind_vars,
      form = form
    )
    if (form == 'ols') {
      fit_model = reg$model
    }
    if (form == 'tobit') {
      fit_model = reg
    }
    wo_fad = test %>%
      mutate(FAD = 0)

    wi_fad = test %>%
      mutate(FAD = 1)

    true_pred_catch =  predict(fit_model,
                               test,
                               res.var = fit_model$VCOV,
                               se.fit = T)$fit


    no_fad_catch = predict(fit_model,
                           wo_fad,
                           res.var = fit_model$VCOV,
                           se.fit = T)$fit

    with_fad_catch = predict(fit_model,
                             wi_fad,
                             res.var = fit_model$VCOV,
                             se.fit = T)$fit

  } #close hurdle if
  if (form == 'hurdle') {
    training$logit_dep_var = as.numeric(training[, dep_var] > min(training[, dep_var]))

    training$norm_dep_var = NA

    training$norm_dep_var = training[,dep_var] #as.numeric(unlist(exp(training[, dep_var])))
    # training$norm_dep_var = as.numeric(unlist(exp(training[, dep_var])))

    # logit_ind_vars = ind_vars[!ind_vars %in% c('factor_country','factor_area')]

    logit_fmla <-
      as.formula(paste("logit_dep_var ~ ", paste(ind_vars, collapse = "+")))

    reg_fmla <-
      as.formula(paste("norm_dep_var ~ ", paste(ind_vars, collapse = "+")))

    gam_reg_fmla = ind_vars

    to_smooth = ind_vars[ind_vars %in% gam_smooth]

    gam_reg_fmla[ind_vars %in% gam_smooth] = paste('s(', to_smooth, ')', sep = '')

    gam_reg_fmla = as.formula(paste(dep_var," ~ ", paste(gam_reg_fmla, collapse = "+")))

    logit_model = glm(logit_fmla ,
                      family = binomial(link = logit),
                      data = training)
        # pos_model_glm = glm(
        #   reg_fmla ,
        #   family = gaussian(link = log),
        #   data = training %>% filter(logit_dep_var == 1)
        # )
    # browser()
    #
    # browser()
    pos_model = gam(
      gam_reg_fmla ,
      family = gaussian(link = log),
      # weights = varIdent(~1 | size_bin),
      data = training %>% filter(logit_dep_var == 1))

    # lmc = lmeControl(niterEM = 5000,msMaxIter = 1000)

    # pos_model2 = mgcv::gamm(
    #   gam_reg_fmla ,
    #   method = 'REML',
    #   control  = lmc,
    #   random = list(factor_year = ~1),
    #   weights = varPower(form =~norm_dep_var),
    #   # family = gaussian(link = log),
    #   data = training %>% filter(logit_dep_var == 1)
    # )

    wo_fad = test %>%
      mutate(FAD = 0)

    wi_fad = test %>%
      mutate(FAD = 1)

    logit_pred_catch =  predict(logit_model,
                                test,
                                type = 'response',
                                se = T)

    pos_pred_catch =  predict(pos_model,
                              test,
                              type = 'response',
                              se = T)
    # test %>%
    #   group_by(country) %>%
    #   summarise(num_fads = length(unique(FAD)))

    #
        # a = augment(logit_model)
        # a$prob = predict(logit_model,
        #                  training,
        #                  type = 'response',
        #                  se = T)$fit

    true_pred_catch = log(logit_pred_catch$fit * pos_pred_catch$fit)

    # (logit_pred_catch$fit * pos_pred_catch$fit)
    #
    # (log(plogis(logit_pred_catch$fit)) + pos_pred_catch$fit)
    # # exp(logit_pred_catch$fit)

    logit_pred_nofad_catch =  predict(logit_model,
                                      wo_fad,
                                      type = 'response',
                                      se = T)

    pos_pred_nofad_catch =  predict(pos_model,
                                    wo_fad,
                                    type = 'response',
                                    se = T)

    no_fad_catch = log(logit_pred_nofad_catch$fit * pos_pred_nofad_catch$fit)

    logit_pred_wifad_catch =  predict(logit_model,
                                      wi_fad,
                                      type = 'response',
                                      se = T)

    pos_pred_wifad_catch =  predict(pos_model,
                                    wi_fad,
                                    type = 'response',
                                    se = T)

    with_fad_catch = log(logit_pred_wifad_catch$fit * pos_pred_wifad_catch$fit)
  }
  if (form == 'constant') {

  }

  preds = data.frame(wi_fad_pred = with_fad_catch,
                     no_fad_pred = no_fad_catch,
                     true_pred = true_pred_catch)

  out <- test %>%
    bind_cols(preds)

  # out_list = list(preds = preds, regs = regs)
  return(out)
}
