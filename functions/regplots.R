#' \code{regplots} is a function for producing standard regression
#' diagnostics from R regression options.
#'
#'
#' @param model a regression object
#' @param fig_dir the figure you want to print plots to
#'
#' @return a list object of regression diagnostics
#' @export
#'
#' @examples regplots(lm(mpg ~ hp, data = mtcars ),fig_dir = 'diagnostics')
#'
regplots <- function(model,fig_dir = ''){

  aug_mod <- augment(model)

# data correlations -------------------------------------------------------

  mod_dat <- model$model %>%
    cor()

# Heteroskedasticity ------------------------------------------------------

resid_vs_fitted_plot = aug_mod %>%
  ggplot(aes(.fitted, .resid)) +
  geom_point(alpha = 0.75) +
  geom_hline(aes(yintercept = 0))


# Distribution ------------------------------------------------------------

resid_hist_plot = aug_mod %>%
  ggplot(aes(.resid)) +
  geom_histogram(bins = 45) +
  geom_vline(aes(xintercept = 0))

resids = rstudent(model)
qqplot(rgamma(length(resids), scale = var(resids, na.rm = T)/mean(resids, na.rm = T),
                      shape = mean(resids, na.rm = T)^2/var(resids, na.rm = T)), resids)
qqnorm(resids)


gg_qq <- function(x, distribution = "norm", ..., line.estimate = NULL, conf = 0.95,
                  labels = names(x)){
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))
  x <- na.omit(x)
  ord <- order(x)
  n <- length(x)
  P <- ppoints(length(x))
  df <- data.frame(ord.x = x[ord], z = q.function(P, ...))

  if(is.null(line.estimate)){
    Q.x <- quantile(df$ord.x, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b <- diff(Q.x)/diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.x ~ z))
  }

  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE

  if(!is.null(labels)){
    df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower, labels[ord],"")
  }

  p <- ggplot(df, aes(x=z, y=ord.x)) +
    geom_point() +
    geom_abline(intercept = coef[1], slope = coef[2]) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2)
  if(!is.null(labels)) p <- p + geom_text( aes(label = label))
  print(p)
  coef
}

gg_qq(x = rstudent(model), distribution = 'gamma')

# Independence (residuals vs. each predictor)

}