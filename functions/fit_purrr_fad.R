#' fit_purrr_fad runs a regression and predicts changes in skj and bet catch from fad reduction
#' using the purrr package
#'
#' @param training data
#' @param test data
#' @param dep_var dependent variable
#' @param ind_vars independent variables
#' @param form type of regression to use ols tobit hurdle
#'
#' @return
#' @export
#'
fit_purrr_fad <- function(training,
                          test,
                          dep_var,
                          ind_vars,
                          form = 'ols',
                          use_gam = T,
                          use_glmer = F,
                          gam_smooth,
                          use_stan = F,
                          num_chains = 1,
                          num_cores = 1,
                          iter = 2000,
                          link_fun = 'identity',
                          reg_fmla = NA,
                          logit_fmla = NA) {
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

    if (is.logical(test$fad) == T) {
      wo_fad = test %>%
        mutate(fad = FALSE)

      wi_fad = test %>%
        mutate(fad = TRUE)

    } else {
      wo_fad <-  test

      wo_fad$fad[wo_fad$fad != 'una'] <- 'una'

      wi_fad <-  test

      wi_fad$fad[wi_fad$fad == 'una'] <- 'afad'

    }

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

    training$norm_dep_var = training[[dep_var]]

    test$logit_dep_var = as.numeric(test[, dep_var] > min(test[, dep_var]))

    test$norm_dep_var = test[[dep_var]]

    if (all(is.na(logit_fmla)) == T)
    {
      logit_fmla <-
        as.formula(paste("logit_dep_var ~ ", paste(ind_vars[!str_detect(ind_vars, '\\|')], collapse = "+")))
    }
    if (all(is.na(reg_fmla)) == T) {
      reg_fmla <-
        as.formula(paste("norm_dep_var ~ ", paste(ind_vars, collapse = "+")))
    }

    gam_reg_fmla = ind_vars

    to_smooth = ind_vars[ind_vars %in% gam_smooth]

    gam_reg_fmla[ind_vars %in% gam_smooth] = paste('s(', to_smooth, ')', sep = '')

    gam_reg_fmla = as.formula(paste(dep_var, " ~ ", paste(gam_reg_fmla, collapse = "+")))
    if (use_stan == F) {
      logit_model = glm(logit_fmla ,
                        family = binomial(link = 'cloglog'),
                        data = training)
    } else {
      logit_model = glm(logit_fmla ,
                        family = binomial(link = 'cloglog'),
                        data = training)
      # logit_model = stan_glm(
      #   logit_fmla ,
      #   family = binomial(link = 'logit'),
      #   data = training,
      #   chains = num_chains,
      #   cores = num_cores
      #)
    }

    if (use_gam == T) {
      pos_model = gam(
        gam_reg_fmla ,
        family = gaussian(link = link_fun),
        data = training %>% filter(logit_dep_var == 1)
      )
    } else {
      if (use_stan == F & use_glmer == F) {
        pos_model = glm(
          reg_fmla ,
          family = gaussian(link = link_fun),
          data = training %>% filter(logit_dep_var == 1)
        )
      }
      if (use_glmer == T) {
        pos_model <-
          lme4::lmer(reg_fmla, data = training %>% filter(logit_dep_var == 1))

      }

      if (use_stan == T) {
        pos_model = stan_glm(
          reg_fmla ,
          family = gaussian(link = 'identity'),
          # family = gaussian(link = link_fun),
          data = training %>% filter(logit_dep_var == 1),
          chains = num_chains,
          cores = num_cores,
          iter = iter,
          adapt_delta = 0.975
        )
      }
    }

    # wo_fad = test %>%
    #       mutate(fad = FALSE)
    #
    #     wi_fad = test %>%
    #       mutate(fad = TRUE)
    #
    if (is.logical(test$fad) == T) {
      wo_fad = test %>%
        mutate(fad = FALSE)

      wi_fad = test %>%
        mutate(fad = TRUE)


    } else {
      wo_fad <-  test

      wo_fad$fad[wo_fad$fad != 'una'] <- 'una'

      wi_fad <-  test

      wi_fad$fad[wi_fad$fad == 'una'] <- 'afad'

    }



    # wi_fad = test %>%
    #   mutate(fad = TRUE)

    shmear <- broom::augment(pos_model) %>%
      ungroup() %>%
      mutate(small_prediction = .fitted <= -2.5) %>%
      group_by(small_prediction) %>%
      summarise(smear_factor = mean(exp(.resid)))

    if (use_stan == F) {
      logit_pred_catch =  predict(logit_model,
                                  test,
                                  type = 'response',
                                  se = T)
      pos_pred_catch =  predict(pos_model,
                                test,
                                type = 'response',
                                se = T, allow.new.levels = T)

      logit_pred_nofad_catch =  predict(logit_model,
                                        wo_fad,
                                        type = 'response',
                                        se = T)
      pos_pred_nofad_catch =  predict(pos_model,
                                      wo_fad,
                                      type = 'response',
                                      se = T,
                                      allow.new.levels = T)


      logit_pred_wifad_catch =  predict(logit_model,
                                        wi_fad,
                                        type = 'response',
                                        se = T)

      pos_pred_wifad_catch =  predict(pos_model,
                                      wi_fad,
                                      type = 'response',
                                      se = T,allow.new.levels = T)


      if (use_glmer == F){

        pos_pred_catch <- pos_pred_catch$fit

        pos_pred_nofad_catch <- pos_pred_nofad_catch$fit

        pos_pred_wifad_catch <- pos_pred_wifad_catch$fit
      }

    } else {
      logit_pred_catch <-   predict(logit_model,
                                    test,
                                    type = 'response',
                                    se = T) %>%
        as_data_frame() #%>%
      # rename(fit = value)

      # pos_pred_catch = predict(pos_model, newdata = test[, c('norm_dep_var', ind_vars)]) %>%
      #   as_data_frame() %>%
      #   rename(fit = value) %>%
      #   mutate(fit = exp(fit))
      pos_pred_catch <-
        posterior_predict(pos_model, newdata = test[, c('norm_dep_var', ind_vars)]) %>%
        as.data.frame() %>%
        gather('entry', 'pred') %>%
        mutate(entry = str_replace(entry, '^.', '') %>% as.numeric()) %>%
        group_by(entry) %>%
        summarise(fit = mean(pred, na.rm = T)) %>%
        ungroup()


      logit_pred_nofad_catch =  predict(logit_model,
                                        wo_fad,
                                        type = 'response',
                                        se = T) %>%
        as_data_frame()

      pos_pred_nofad_catch <-
        posterior_predict(pos_model, newdata = wo_fad[, c('norm_dep_var', ind_vars)]) %>%
        as.data.frame() %>%
        gather('entry', 'pred') %>%
        mutate(entry = str_replace(entry, '^.', '') %>% as.numeric()) %>%
        group_by(entry) %>%
        summarise(fit = mean(pred, na.rm = T)) %>%
        ungroup()


      logit_pred_wifad_catch =  predict(logit_model,
                                        wi_fad,
                                        type = 'response',
                                        se = T) %>%
        as_data_frame()

      pos_pred_wifad_catch =  posterior_predict(pos_model, newdata = wi_fad[, c('norm_dep_var', ind_vars)]) %>%
        as.data.frame() %>%
        gather('entry', 'pred') %>%
        mutate(entry = str_replace(entry, '^.', '') %>% as.numeric()) %>%
        group_by(entry) %>%
        summarise(fit = mean(pred, na.rm = T)) %>%
        ungroup()

    }

    # Calculate smearing estimate
    variance_check <- pos_model %>%
      augment() %>%
      ggplot(aes(.fitted, .resid)) +
      geom_point()

    shmear <- augment(pos_model) %>%
      ungroup() %>%
      mutate(small_prediction = .fitted <= -2.5) %>%
      group_by(factor_area, fad) %>%
      summarise(smear_factor = mean(exp(.resid)))

      pos_predictions <-
      data_frame(
        factor_area = test$factor_area,
        fad = test$fad,
        pos_pred_catch = as.numeric(pos_pred_catch),
        pos_pred_nofad_catch = as.numeric(pos_pred_nofad_catch),
        pos_pred_wifad_catch = as.numeric(pos_pred_wifad_catch)
      ) %>%
      mutate(obs = 1:dim(.)[1]) %>%
      gather('blah', 'fit', contains('pos_pred')) %>%
      mutate(small_prediction = fit <= -2.5) %>%
      left_join(shmear, by = c('factor_area', 'fad')) %>%
      # left_join(shmear, by = 'small_prediction') %>%
      mutate(unbiased_fit = exp(fit) * smear_factor) %>%
      select(obs, blah, unbiased_fit) %>%
      ungroup() %>%
      spread(blah, unbiased_fit) %>%
      select(-obs)

    logit_predictions <-
      data_frame(
        logit_pred_catch = as.numeric(logit_pred_catch$fit),
        logit_pred_nofad_catch = as.numeric(logit_pred_nofad_catch$fit),
        logit_pred_wifad_catch = as.numeric(logit_pred_wifad_catch$fit)
      )


    predictions <- (pos_predictions * logit_predictions) %>%
      as_data_frame()

    true_pred_catch <-  predictions$pos_pred_catch

    no_fad_catch <-  predictions$pos_pred_nofad_catch

    with_fad_catch <-  predictions$pos_pred_wifad_catch

  }
  if (form == 'constant') {

  }
  preds = data.frame(wi_fad_pred = with_fad_catch,
                     no_fad_pred = no_fad_catch,
                     true_pred = true_pred_catch)

  out <- test %>%
    bind_cols(preds)

  return(list(
    pos_model = pos_model,
    logit_model = logit_model,
    predictions = out
  ))
}
