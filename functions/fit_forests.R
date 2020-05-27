#' fit_forests fits ensembles of random forests and other
#' candidate models
#'
#' @param dat
#' @param reg_year
#' @param tree_vars
#' @param ind_vars
#' @param iterations
#' @param glmer_fmla
#' @param dep_var
#' @param form
#' @param run_in_parallel
#' @param num_cores
#' @param use_hurdle
#'
#' @return a list of fitted models
#' @export
#'
fit_forests <- function(dat,
                        reg_year,
                        test_year = 2009,
                        tree_vars,
                        ind_vars,
                        iterations,
                        glmer_fmla,
                        dep_var,
                        form,
                        run_in_parallel = F,
                        num_cores = 1,
                        use_hurdle = T,
                        trees = 2000,
                        holdout_mode = 1) {
  
  forest_dat <- dat %>%
    filter(year >= reg_year,!(
      factor_country %in% c('China', 'Japan', 'South Korea', 'New Zealand')
    )) %>%
    mutate(any_seen = (log_cpue > min(log_cpue, na.rm = T)) %>% as.factor())

  
  used_rows <- forest_dat$rownames

  if (use_hurdle == T){
    seen_dat <- forest_dat %>%
      filter(any_seen == T)

  } else {
  seen_dat <- forest_dat
  }

  if (run_in_parallel == T) {
    doMC::registerDoMC(cores = num_cores)

  }
# 
#   tree_vars <-  c('fad',
#                   'factor_area',
#                   'factor_country',
#                   'sst',
#                   'month',
#                   'year')
#   
  if (holdout_mode == 1){
  
  tree_training <- forest_dat[!(forest_dat$year >= test_year & forest_dat$month %in% c(7,8,9)), c('cpue',tree_vars)]
  
  tree_testing <- forest_dat[(forest_dat$year >= test_year & forest_dat$month %in% c(7,8,9)), c('cpue',tree_vars)]
  
  } else {
    
    tree_training <- forest_dat[(forest_dat$year < test_year), c('cpue',tree_vars)]
    
    tree_testing <- forest_dat[(forest_dat$year >= test_year), c('cpue',tree_vars)]
  }
  
  # tree_training <- forest_dat[forest_dat$year < test_year, c('cpue',tree_vars)]
  # 
  # tree_testing <- forest_dat[forest_dat$year >= test_year, c('cpue',tree_vars)]
  # 
  
  # training <- forest_dat[forest_dat$year < test_year, ]
  # 
  # testing <- forest_dat[forest_dat$year >= test_year, ]
  # 
  # 
  tree_training$date <- lubridate::ymd(paste(tree_training$year,tree_training$month,"01"))
  
  # tree_training <- tree_training %>% 
  #   select(-year)
  # 
  forest_rolling_origin <- rsample::rolling_origin(tree_training %>% nest(-date),
                                                   initial = 90,
                                                   assess = 12,
                                                   cumulative = TRUE)

  tree_training <- tree_training %>% 
    select(-date)
  # tune model
  
  tune_grid <- cross_df(list(
    splitrule = c("extratrees"),
    mtry = seq(2,8, by = 2),
    min_n = c(2,5),
    id = unique(forest_rolling_origin$id)
  )) %>%
    left_join(forest_rolling_origin, by = "id")
  
  future::plan(future::multiprocess, workers = num_cores)
  
  tune_grid <- tune_grid %>%
    mutate(tuning_fit = furrr::future_pmap(
      list(
        mtry = mtry,
        splitrule = splitrule,
        min_n = min_n,
        split = splits
      ),
      tune_forest,
      .progress = TRUE,
      model_class = "rand_forest",
      trees = trees
    ))
  
  best_params <- tune_grid %>% 
    select(-splits) %>% 
    unnest() %>% 
    group_by(mtry, splitrule, min_n) %>% 
    yardstick::rmse(observed, .pred) 
  
  best_params %>%
    ggplot(aes(mtry, .estimate, color = splitrule)) +
    geom_point() +
    facet_wrap(~min_n)

  best <- best_params %>% 
    filter(.estimate == min(.estimate))
  
  forest_recipe <- recipe(cpue ~ ., data = tree_training) %>% 
    step_dummy(all_nominal()) %>% 
    step_nzv(-all_outcomes()) %>% 
    step_corr(-all_outcomes()) %>% 
    prep(data = tree_training, retain = TRUE)
  
  trained_forest <-
    parsnip::rand_forest(
      mode = "regression",
      mtry = best$mtry,
      min_n = best$min_n,
      trees = trees
    ) %>%
    parsnip::set_engine(
      "ranger",
      importance = "none",
      splitrule = best$splitrule
    ) %>%
    parsnip::fit(formula(forest_recipe), data = juice(forest_recipe))
  
  # trained_tree <-
  #   parsnip::boost_tree(
  #     mode = "regression",
  #     mtry = best$mtry,
  #     min_n = best$min_n,
  #     trees = 1000
  #   ) %>%
  #   parsnip::set_engine(
  #     "xgboost") %>%
  #   parsnip::fit(formula(forest_recipe), data = juice(forest_recipe))
  
  pred <- predict(trained_forest, new_data = bake(forest_recipe,forest_dat))
  
  # pred <- predict(trained_tree, new_data = bake(forest_recipe,forest_dat))
  
  forest_predictions <- forest_dat %>%
    mutate(
      log_cpue_hat =  pred$.pred, #ifelse(use_hurdle == T,pos_pred, pos_pred),
      catch_hat =  pmax(0,(log_cpue_hat)) * days,
      cpue_hat =   pmax(0,(log_cpue_hat))
    ) 

# browser()
#   forest_predictions %>%
#     mutate(testing = year >= 2009) %>%
#     ggplot(aes(cpue, cpue_hat, color = year)) +
#     geom_point() +
#     geom_smooth(method = "lm") +
#     geom_abline(aes(intercept = 0, slope = 1))+
#     facet_wrap(~testing)

  
  # train_test_frame <- tibble(training = list(training), test = list(testing))
  
  # glm_model <-train_test_frame %>% 
  #   # pblapply(
  #   #   'all',
  #   #   partition_tuna,
  #   #   dat = dat %>% filter(year >= reg_year),
  #   #   test_set = 1:dim(dat %>% filter(year >= reg_year))[1]
  #   # ) %>%
  #   # purrr::transpose() %>% #flatten data
  #   # as_data_frame() %>% #convert to tibble
  #   mutate(
  #     fitted_model = map2(
  #       training,
  #       test,
  #       fit_purrr_fad,
  #       iter = iterations,
  #       reg_fmla = NA,
  #       logit_fmla = NA,
  #       dep_var = dep_var,
  #       ind_vars = ind_vars,
  #       form = form,
  #       gam_smooth = F,
  #       use_gam = F,
  #       use_glmer = F,
  #       use_stan = F,
  #       num_chains = 1,
  #       num_cores = 1,
  #       link_fun = 'identity'
  #     )
  #   )

  # glmer_model <- train_test_frame %>%
  #   mutate(
  #     fitted_model = map2(
  #       training,
  #       test,
  #       fit_purrr_fad,
  #       iter = iterations,
  #       reg_fmla = glmer_fmla,
  #       logit_fmla = NA,
  #       dep_var = dep_var,
  #       ind_vars = ind_vars,
  #       form = form,
  #       gam_smooth = F,
  #       use_gam = F,
  #       use_glmer = T,
  #       use_stan = F,
  #       num_chains = 1,
  #       num_cores = 1,
  #       link_fun = 'identity'
  #     )
  #   )

  # glm_predictions <-
  #   exp(predict(glm_model$fitted_model[[1]]$pos_model, newdata = forest_dat)) *
  #   
  #   predict(glm_model$fitted_model[[1]]$logit_model,
  #           newdata = forest_dat,
  #           type = "response")
  # 
  # glm_predictions <- forest_dat %>% 
  #   mutate(true_pred = glm_predictions) %>% 
  #   select(year, month, catch, true_pred, days, cpue) %>%
  #   mutate(catch_hat = true_pred * days,
  #          method = 'glm')

  
  # glmer_predictions <- predict(glmer_model$fitted_model[[1]]$pos_model, newdata = forest_dat) * 
  #   
  #   predict(glmer_model$fitted_model[[1]]$logit_model, newdata = forest_dat, type = "response")
  # 
  # 
  # glmer_predictions <- glmer_model$fitted_model[[1]]$predictions %>%
  #   select(year, month, catch, true_pred, days) %>%
  #   mutate(catch_hat = true_pred * days,
  #          method = 'glmer')
  # 
  # glm_predictions <- glm_model$fitted_model[[1]]$predictions %>%
  #   select(year, month, catch, true_pred, days) %>%
  #   mutate(catch_hat = true_pred * days,
  #          method = 'glm')
  # 
  # glmer_predictions <- glmer_model$fitted_model[[1]]$predictions %>%
  #   select(year, month, catch, true_pred, days) %>%
  #   mutate(catch_hat = true_pred * days,
  #          method = 'glmer')

  forest_preds <- forest_predictions %>%
    select(year, month, catch, log_cpue_hat, days, catch_hat, cpue) %>%
    dplyr::rename(true_pred = log_cpue_hat) %>%
    mutate(method = 'random-forest')
  
  model_comp <- forest_preds #%>%
    # bind_rows(glmer_predictions) %>%
    # bind_rows(glm_predictions)
  
  # model_comp %>% 
  #   mutate(testing = year >= 2009) %>% 
  #   ggplot(aes(cpue, true_pred)) + 
  #   geom_point() + 
  #   geom_smooth(method = "lm") + 
  #   geom_abline(aes(intercept = 0, slope = 1))+
  #   facet_grid(method~testing)
  # 
  

  model_comp <- model_comp %>%
    mutate(
      residuals = catch_hat - catch,
      squared_error = residuals ^ 2,
      year_month = zoo::as.yearmon(paste(year, month, sep = '-'))
    )

  model_comp <- model_comp %>%
    gather('catch_source', 'catch', catch, catch_hat)


  return(
    list(
      model_comp = model_comp,
      # logit_forest = logit_fit,
      seen_forest = trained_forest,
      forest_predictions = forest_predictions,
      forest_recipe = forest_recipe
    )
  )

}
