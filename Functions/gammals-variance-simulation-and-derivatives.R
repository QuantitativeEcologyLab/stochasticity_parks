#' this code was originally written by Gavin Simpson. You may want
#' to look at the `{gratia}` package that Gavin maintains, since he
#' intends to include all of these functions in the package. The functions
#' for derivatives for generic location (H)GAMs aleady exist, but he was
#' working on writing generic functions for location-scale families, too.
#' *original code*: https://github.com/simpson-lab/wpg-mb-lakes/blob/main/analysis/variance-simulation-and-derivatives.R
EXAMPLES <- FALSE # should examples be run?

##' Simulate from the posterior distribution of the mean of a Gamma LS GAM
##'  using Gaussian approximation to the posterior
##'
##' @param model         the fitted GAM
##' @param data          the new data locations you want to get variance for
##' @param nsims         the number of posterior draws wanted - low by default
##'                      to avoid excessive computation, but needs to be
##'                      10,000+ for quantile-based intervals
##' @param unconditional logical; use the smoothness selection corrected version
##'                      of the Bayesian covariance matrix of the model?
`gammals_mean` <- function(model, data, nsims = 100, terms = NULL,
                           unconditional  = FALSE, ...) {
  ## Simulate variance from posterior
  sim <- sim_gammals_mean(model = model, data = data, terms = terms,
                          nsims = nsims, unconditional = unconditional)
  ## process results into a tibble
  colnames(sim) <- paste0("sim", seq_len(nsims))
  tbl <- as_tibble(sim) %>%
    bind_cols(data) %>%
    tibble::rowid_to_column(var = "row")
  tbl <- pivot_longer(tbl,
                      cols = matches("^sim"),
                      names_to = "simulation",
                      values_to = "mean")
  tbl
}

##' Simulate from the posterior distribution of the variance of a Gamma LS GAM
##'  using Gaussian approximation to the posterior
##'
##' @param model         the fitted GAM
##' @param data          the new data locations you want to get variance for
##' @param nsims         the number of posterior draws wanted - low by default
##'                      to avoid excessive computation, but needs to be
##'                      10,000+ for quantile-based intervals
##' @param unconditional logical; use the smoothness selection corrected version
##'                      of the Bayesian covariance matrix of the model?
`gammals_var` <- function(model, data, nsims = 100, terms = NULL,
                          unconditional  = FALSE, ...) {
  ## Simulate variance from posterior
  sim <- sim_gammals_var(model = model, data = data, terms = terms,
                         nsims = nsims, unconditional = unconditional)
  ## process results into a tibble
  colnames(sim) <- paste0("sim", seq_len(nsims))
  tbl <- as_tibble(sim) %>%
    bind_cols(data) %>%
    tibble::rowid_to_column(var = "row")
  tbl <- pivot_longer(tbl,
                      cols = matches("^sim"),
                      names_to = "simulation",
                      values_to = "variance")
  tbl
}

##' Simulate from the posterior distribution of the derivative of the mean
##' of a location -scale Gamma GAM using Gaussian approximation to the
##' posterior. Derivatives are computed using left finite differences only.
##' Other arguments as per above.
##'
##' @param var character; the variable to shift along for derivatives
##' @param eps numeric; the value to shift data by for finite differences
`gammals_mean_deriv` <- function(model, data, var, nsims = 100, terms = NULL,
                                 unconditional  = FALSE, eps = 1e-07, ...) {
  
  ## f'(x) = (f(x + eps) - f(x))/eps as eps --> 0
  
  ## prediction matrix
  Xp1 <- predict(model, newdata = data, type = 'lpmatrix', terms = terms)
  ## model parameters
  coefs <- coef(model)
  ## Bayesian covariance matrix
  Vb <- vcov(model, unconditional = unconditional)
  ## which coefs go with the theta linear predictor
  theta_take <- grepl('^s\\.1', colnames(Xp1)) |
    colnames(Xp1) %in% c('(Intercept).1')
  mu_take <- !theta_take
  ## inverse link functions
  ilink_mu <- inv_link(model, parameter = "location") # mu inv link function
  
  ## Simulate from posterior using Gaussian approximation
  betas <- mvnfast::rmvn(n = nsims,
                         mu = coefs,
                         sigma = Vb)
  
  ## predict variance for data
  mu1 <- est_gammals_mean(betas, Xp1, mu_take, ilink_mu)
  data2 <- data # copy
  ## shift the variable of interest by eps
  data2[[var]] <- data2[[var]] + eps
  ## predict for shifted data
  ## prediction matrix
  Xp2 <- predict(model, newdata = data2, type = 'lpmatrix', terms = terms)
  mu2 <- est_gammals_mean(betas, Xp2, mu_take, ilink_mu)
  
  ## compute finite differences
  sim_d <- (mu2 - mu1) / eps
  
  ## process into a tibble
  colnames(sim_d) <- paste0("sim", seq_len(nsims))
  tbl <- as_tibble(sim_d) %>%
    bind_cols(data) %>%
    tibble::rowid_to_column(var = "row")
  tbl <- pivot_longer(tbl, cols = matches("^sim"),
                      names_to = "simulation",
                      values_to = "derivative")
  tbl
}

##' Simulate from the posterior distribution of the derivative of the variance
##' of a location -scale Gamma GAM using Gaussian approximation to the
##' posterior. Derivatives are computed using left finite differences only.
##' Other arguments as per above.
##'
##' @param var character; the variable to shift along for derivatives
##' @param eps numeric; the value to shift data by for finite differences
`gammals_var_deriv` <- function(model, data, var, nsims = 100, terms = NULL,
                                unconditional  = FALSE, eps = 1e-07, ...) {
  
  ## f'(x) = (f(x + eps) - f(x))/eps as eps --> 0
  
  ## prediction matrix
  Xp1 <- predict(model, newdata = data, type = 'lpmatrix', terms = terms)
  ## model parameters
  coefs <- coef(model)
  ## Bayesian covariance matrix
  Vb <- vcov(model, unconditional = unconditional)
  ## which coefs go with the theta linear predictor
  theta_take <- grepl('^s\\.1', colnames(Xp1)) |
    colnames(Xp1) %in% c('(Intercept).1')
  mu_take <- !theta_take
  ## inverse link functions
  ilink_mu <- inv_link(model, parameter = "location") # mu inv link function
  ilink_theta <- inv_link(model, parameter = "scale") # theta inverse link
  
  ## Simulate from posterior using Gaussian approximation
  betas <- mvnfast::rmvn(n = nsims,
                         mu = coefs,
                         sigma = Vb)
  
  ## predict variance for data
  var1 <- est_gammals_var(betas, Xp1, mu_take, theta_take,
                          ilink_mu, ilink_theta)
  data2 <- data # copy
  ## shift the variable of interest by eps
  data2[[var]] <- data2[[var]] + eps
  ## predict for shifted data
  ## prediction matrix
  Xp2 <- predict(model, newdata = data2, type = 'lpmatrix', terms = terms)
  var2 <- est_gammals_var(betas, Xp2, mu_take, theta_take,
                          ilink_mu, ilink_theta)
  
  ## compute finite differences
  sim_d <- (var2 - var1) / eps
  
  ## process into a tibble
  colnames(sim_d) <- paste0("sim", seq_len(nsims))
  tbl <- as_tibble(sim_d) %>%
    bind_cols(data) %>%
    tibble::rowid_to_column(var = "row")
  tbl <- pivot_longer(tbl, cols = matches("^sim"),
                      names_to = "simulation",
                      values_to = "derivative")
  tbl
}

`est_gammals_mean` <- function(betas, Xp, mu_take, ilink_mu) {
  ## subset Xp matrix into mean part
  Xp_mu <- Xp[, mu_take, drop = FALSE]
  
  ## Predict for mean
  fit_mu <- Xp_mu %*% t(betas[, mu_take, drop = FALSE]) # predict on internal scale
  fit_mu <- ilink_mu(fit_mu) # apply g-1() this is just identity so redundant
  ## This model is parameterised in terms of the log-mean, so we still need to
  ## transform to the actual data scale using exp()
  fit_mu <- exp(fit_mu)
  fit_mu
}


`est_gammals_var` <- function(betas, Xp, mu_take, theta_take,
                              ilink_mu, ilink_theta) {
  ## subset Xp matrix into mean and scale parts
  Xp_mu <- Xp[, mu_take, drop = FALSE]
  Xp_theta <- Xp[, theta_take, drop = FALSE]
  
  ## Predict for mean
  fit_mu <- Xp_mu %*% t(betas[, mu_take, drop = FALSE]) # predict on internal scale
  fit_mu <- ilink_mu(fit_mu) # apply g-1() this is just identity so redundant
  ## This model is parameterised in terms of the log-mean, so we still need to
  ## transform to the actual data scale using exp()
  fit_mu <- exp(fit_mu)
  
  ## Predict for theta
  fit_theta <- Xp_theta %*%
    t(betas[, theta_take, drop = FALSE]) # predict on internal scale
  fit_theta <- ilink_theta(fit_theta) # apply g-1()
  ## fit theta even after using inverse link is log for theta parameter, so
  ## more back transforming
  fit_theta <- exp(fit_theta)
  
  ## variance is mu * s where s is the scale in sense of rgamma, not theta
  ## From ?rgamma Var(y) = shape * scale^2 = (1 / theta) * (mu * theta)^2
  ## From ?gammals Var(y) = mu * scale = mu * s
  ##   where scale = s = mu * theta. Hence from ?gammals we arrive finally
  ##   at: Var(y) = mu * s = mu * (mu * theta)
  fit_var_draws <- fit_mu * (fit_mu * fit_theta)
  ## return
  fit_var_draws
}

##' The internal workhorse does all the cool stuff 
`sim_gammals_mean` <- function(model, data, nsims = 100, terms = NULL,
                               unconditional = FALSE, ...) {
  ## prediction matrix
  Xp <- predict(model, newdata = data, type = 'lpmatrix', terms = terms)
  ## model parameters
  coefs <- coef(model)
  ## Bayesian covariance matrix
  Vb <- vcov(model, unconditional = unconditional)
  ## which coefs go with the theta linear predictor
  theta_take <- grepl('^s\\.1', colnames(Xp)) |
    colnames(Xp) %in% c('(Intercept).1')
  
  ## Simulate from posterior using Gaussian approximation
  betas <- mvnfast::rmvn(n = nsims,
                         mu = coefs,
                         sigma = Vb)
  
  ## simplify later code so form the compliment to select mean
  ## linear predictor
  mu_take <- !theta_take
  ## subset Xp matrix into mean and theta parts
  Xp_mu <- Xp[, mu_take, drop = FALSE]
  
  ## Predict for mean
  fit_mu <- Xp_mu %*% t(betas[, mu_take, drop = FALSE]) # predict on internal scale
  ilink_mu <- inv_link(model, parameter = "location") # link function
  fit_mu <- ilink_mu(fit_mu) # apply g-1() this is just identity so redundant
  ## This model is parameterised in terms of the log-mean, so we still need to
  ## transform to the actual data scale using exp()
  fit_mu <- exp(fit_mu)
  fit_mu
}

##' The internal workhorse does all the cool stuff 
`sim_gammals_var` <- function(model, data, nsims = 100, terms = NULL,
                              unconditional = FALSE, ...) {
  ## prediction matrix
  Xp <- predict(model, newdata = data, type = 'lpmatrix', terms = terms)
  ## model parameters
  coefs <- coef(model)
  ## Bayesian covariance matrix
  Vb <- vcov(model, unconditional = unconditional)
  ## which coefs go with the theta linear predictor
  theta_take <- grepl('^s\\.1', colnames(Xp)) |
    colnames(Xp) %in% c('(Intercept).1')
  
  ## Simulate from posterior using Gaussian approximation
  betas <- mvnfast::rmvn(n = nsims,
                         mu = coefs,
                         sigma = Vb)
  
  ## simplify later code so form the compliment to select mean
  ## linear predictor
  mu_take <- !theta_take
  ## subset Xp matrix into mean and theta parts
  Xp_mu <- Xp[, mu_take, drop = FALSE]
  Xp_theta <- Xp[, theta_take, drop = FALSE]
  
  ## Predict for mean
  fit_mu <- Xp_mu %*% t(betas[, mu_take, drop = FALSE]) # predict on internal scale
  ilink_mu <- inv_link(model, parameter = "location") # link function
  fit_mu <- ilink_mu(fit_mu) # apply g-1() this is just identity so redundant
  ## This model is parameterised in terms of the log-mean, so we still need to
  ## transform to the actual data scale using exp()
  fit_mu <- exp(fit_mu)
  
  ## Predict for theta
  fit_theta <- Xp_theta %*%
    t(betas[, theta_take, drop = FALSE]) # predict on internal scale
  ilink_theta <- inv_link(model, parameter = "scale") # theta inverse link
  fit_theta <- ilink_theta(fit_theta) # apply g-1()
  ## fit scale even after using inverse link is log for theta parameter, so
  ## more back transforming
  fit_theta <- exp(fit_theta)
  
  ## variance is mu * s where s is the scale in sense of rgamma, not theta
  ## From ?rgamma Var(y) = shape * scale^2 = (1 / theta) * (mu * theta)^2
  ## From ?gammals Var(y) = mu * scale = mu * s
  ##   where scale = s = mu * theta. Hence from ?gammals we arrive finally
  ##   at: Var(y) = mu * s = mu * (mu * theta)
  fit_var_draws <- fit_mu * (fit_mu * fit_theta)
  ## return
  fit_var_draws
}