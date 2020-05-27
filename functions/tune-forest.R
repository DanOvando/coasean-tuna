#' tune forest
#'
#' @param split the data split
#' @param mtry mtry
#' @param splitrule splitrule
#' @param min_n minimum node size
#'
#' @return observed and predicted
#' @export
#'
tune_forest <- function(split, mtry, splitrule, min_n,
                        trees = 2000,
                        model_class = "rand_forest"
) {
  
  analysis_split <-   rsample::analysis(split) %>%
    unnest() %>% 
    select(-date)
  
  assessment_split <-  rsample::assessment(split) %>%
    unnest() %>% 
    select(-date)
  
  forest_recipe <- recipe(cpue ~ ., data = analysis_split) %>% 
    step_dummy(all_nominal()) %>% 
    step_nzv(-all_outcomes()) %>% 
    step_corr(-all_outcomes()) %>% 
    prep(data = split, retain = TRUE)
  
  if (model_class == "rand_forest"){
  
  forest_assess <-
    parsnip::rand_forest(
      mode = "regression",
      mtry = mtry,
      min_n = min_n,
      trees = trees
    ) %>%
    parsnip::set_engine(
      "ranger",
      importance = "none",
      splitrule = splitrule) %>%
    parsnip::fit(formula(forest_recipe), data = juice(forest_recipe))
  } else if (model_class == "boost_tree") {
    
    forest_assess <-
      parsnip::boost_tree(
        mode = "regression",
        mtry = mtry,
        min_n = min_n,
        trees = trees
      ) %>%
      parsnip::set_engine(
        "xgboost",
        importance = "none") %>%
      parsnip::fit(formula(forest_recipe), data = juice(forest_recipe))
    
    
  }
  
  # importance(forest_assess$fit) %>% 
  #   broom::tidy() %>% 
  #   View()
  
  assessment_prediction <-
    predict(forest_assess, new_data = bake(forest_recipe, assessment_split))
  
  assessment_prediction$observed = assessment_split$cpue
  
  return(assessment_prediction)
  
} # close tune forest
