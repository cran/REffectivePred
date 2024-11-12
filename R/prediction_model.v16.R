#' @importFrom graphics abline axis lines points rect
#' @importFrom stats dlnorm dnorm optim
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics title
NULL

#' @title Detect start of waves
#' @description Find the approximate beginning times of waves given times of peaks. This is based on local maxima of \eqn{R_t}.
#'
#' @param rt_values A vector containing rt values for each time point.
#' @param peaks_x Time points of peaks based on daily cases (one for each wave).
#' @param search_range The range of data points to go through for filtering invalid beginnings (a vector).
#' @details Note: This is provided for convenience only, and is not meant to replace an analyst's determination of wave bounds.
#' @returns A list containing the detected start times of waves, as scalars.
#'
#' @export

find_starts <- function(rt_values, peaks_x, search_range) {
  max_positions <- c()

  # go through all given peaks
  for (i in 1:length(peaks_x)) {
    peak_x <- peaks_x[i]
    max_position <- peak_x - 1

    actual_search_range <- search_range

    if ((peak_x <= search_range) || (peak_x - max_positions[length(max_positions)] <= search_range)) {
      if (i == 1) {
        # if first peak too early
        actual_search_range <- peak_x - 1
      } else {
        # if the distance between two peaks is smaller than the desired search range, divide it by 2
        actual_search_range <- ((peak_x - peaks_x[i - 1]) / 2) - 1
      }
    } else {
      actual_search_range <- search_range
    }

    for (j in (peak_x - 1):(peak_x - actual_search_range)) {
      if (!is.na(rt_values[j]) && rt_values[j] > rt_values[max_position]) {
        max_position <- j
      }
    }
    max_positions <- c(max_positions, max_position)
  }

  output <- max_positions
}


#' @title Detect end of waves
#' @description Find the approximate end times of waves given times of peaks. This is based on local minima of \eqn{R_t}.
#'
#' @param rt_values A vector containing rt values for each time point.
#' @param peaks_x Time points of peaks based on daily cases (one for each wave).
#' @param search_range The range of data points to go through for filtering invalid ends (a vector).
#' @details Note: This is provided for convenience only, and is not meant to replace an analyst's determination of wave bounds.
#' @returns A list containing the detected end times of waves, as scalars.
#'
#' @export

find_ends <- function(rt_values, peaks_x, search_range) {
  min_positions <- c()

  for (i in 1:length(peaks_x)) {
    peak_x <- peaks_x[i]
    min_position <- peak_x + 1

    actual_search_range <- search_range

    # if the distance between two peaks is smaller than the desired search range, divide it by 2
    if ((i != length(peaks_x)) && ((peaks_x[i + 1] - peaks_x[i]) <= search_range)) {
      actual_search_range <- ((peaks_x[i + 1] - peaks_x[i]) / 2) - 1
    } else if ((i == length(peaks_x)) && (peak_x >= (length(rt_values) - search_range))) {
      actual_search_range <- length(rt_values) - peak_x
    }

    for (j in (peak_x + 1):(peak_x + actual_search_range)) {
      if (!is.na(rt_values[j]) && rt_values[j] < rt_values[min_position]) {
        min_position <- j
      }
    }
    min_positions <- c(min_positions, min_position)
  }

  output <- min_positions
}


#' @title Utility function for range manipulation
#' @description Converts a list of waves to a two-dimensional list.
#'
#' @param waves_list A list containing ranges (start, end) of each wave.
#' @returns A two-dimensional list with individual waves as sub-lists.
#' @export

ranges_to_waves <- function(waves_list) {
  tmp_waves_list <- list()
  for (i in c(1:length(waves_list))) {
    tmp_waves_list <- append(tmp_waves_list, list(c(waves_list[[i]][[1]]:waves_list[[i]][[2]])))
  }
  waves_list <- tmp_waves_list
  output <- waves_list
}


#' @title Utility function for range manipulation
#' @description Combines multiple waves into one vector.
#'
#' @param num_waves Total number of waves.
#' @param waves_list A list containing individual waves as sub-lists.
#' @returns A vector with the first \code{num_waves} waves combined.
#' @export


waves_1d_list <- function(num_waves, waves_list) {
  waves <- c()
  for (i in 1:num_waves) {
    waves <- c(waves, unlist(waves_list[i]))
  }
  # waves_pred <- c(waves, pred_save)

  output <- waves
}


#' @title Empirical estimate of \eqn{R_t}
#' @description Compute empirical \eqn{R_t}, via Cori et al. (2013) method.
#'
#' @param cases Vector of (confirmed) cases.
#' @param window_size The maximum value for the serial interval.
#' @param serial_mean Mean of the serial interval on the log scale.
#' @param serial_var Variance of the serial interval on the log scale.
#' @returns A vector of same length as cases, giving the empirical estimate of the effective reproductive number over time.
#' @export

rt_empirical <- function(cases, window_size, serial_mean, serial_var) {
  cases_extended <- c(rep(0, window_size), cases)
  lambda_extended <- 0 * cases_extended
  for (j in (window_size + 1):length(cases_extended)) {
    lambda_extended[j] <- sum(cases_extended[(j - window_size):(j - 1)] * serial.helper(window_size = window_size, serial_mean = serial_mean, serial_var = serial_var))
  }
  rt_extended <- cases_extended / lambda_extended
  rt.empirical <- rt_extended[(window_size + 1):length(rt_extended)]
  rt.empirical[1] <- NA

  output <- rt.empirical
}


#' @title The likelihood function
#' @description The negative log likelihood function of the model.
#'
#' @param param Includes the following sets of parameters in a vector, in this order:
#'  \itemize{
#'    \item a1,a2,a3,a4 - Parameters for curve c() specified by \code{rt_func}.
#'    \item nu - Loss of Immunity rate.
#'    \item v2,v3,v4,v5 - Transmissibility of variants in waves 2+, as relative multiplication factors compared to transmissibility in wave 1.
#'    \item psi1,psi2,psi3,psi4 - Psi parameters for severity levels 1,2,3 and 4.
#'    \item u,v - Variance parameters. Only u is currently in use.
#'    \item beta0,beta.R,beta.E,beta.W - When restrictions = NULL. Currently unsupported.
#'  }
#' @param params_limits Boundaries/limits of the ini_params.
#' @param restrictions A numeric integer vector giving the severity of restrictions.
#'   Zero means no restriction, and higher numbers means greater severity/disruption.
#'   The ordered unique values should be consecutive integers starting from zero.
#'   Each number (other than 0) adds a new parameter to the fit. restrictions = NULL
#'   causes the function to use mobility data instead of the psi values (currently unsupported).
#' @param restriction.starts A vector of same length as restrictions, of times when restrictions
#'   came into effect. Note: the first index time should be 1.
#' @param ranges An vector of time ranges for the different waves.
#'   The wave ranges should be contiguous, with at least one unit of time
#'   between consecutive waves.
#' @param rt_func The parametric form of function c(). Options are listed under function c_helper.
#' @param silence.errors  Ignores (skips) NA or NaN values when summing up likelihood contributions over time.
#' @param fit.t.pred Time of prediction.
#' @param lt Length of cases.
#' @param cases A vector containing cases for each time-point.
#' @param scenario A character string describing options to deal with restrictions. Currently unsupported.
#' @param H.E Mobility metrics for category Retail & Entertainment. Currently unsupported.
#' @param H.W Mobility metrics for category Workplaces. Currently unsupported.
#' @param adj.period Delays in society adjusting.
#' @param population total population size.
#' @param rho Under-reporting fraction.
#' @param serial_mean Mean of the serial interval on the log scale.
#' @param serial_var Variance of the serial interval on the log scale.
#' @param window_size The maximum value for the serial interval.
#' @returns The negative log likelihood value of the data.
#' @details The predicted curve is computed based on parameters supplied, by first calling the prediction
#' function \code{pred.curve}. The probability model used to compute the likelihood assumes that observed
#' infection at time \eqn{t} are \eqn{\sim N(mean = I_t, sd = \sqrt{u*I_t})}, where \eqn{I_t} are predicted infections, and sums the
#' log-likelihood contributions for each time \eqn{t} during waves, and up to \code{fit.t.pred}.
#' @export

log_lklh <- function(param,
                     params_limits,
                     restrictions = NULL,
                     restriction.starts = NULL,
                     ranges = NULL,
                     rt_func = 1,
                     silence.errors = FALSE,
                     fit.t.pred,
                     lt,
                     cases,
                     scenario = NULL,
                     H.E,
                     H.W,
                     adj.period,
                     population,
                     rho,
                     serial_mean,
                     serial_var,
                     window_size) {
  infections <- cases / rho # Estimated time series of total infections
  yt <- infections

  param <- modif.helper(param, params_limits)

  # Note: indexing from 1 is required  for !is.null(restrictions).

  yt.pr <- NA * yt
  rt.pr <- NA * yt
  lambdat.pr <- NA * yt
  yt.pred <- pred.curve(
    a1 = param[1],
    a2 = param[2],
    a3 = param[3],
    a4 = param[4],
    nu = param[5],
    Psi = param[10:15],
    betas = param[16:19],
    restrictions = restrictions,
    restriction.starts = restriction.starts,
    ranges = ranges,
    variant.transm = c(1, param[6:9]),
    rt_func = rt_func,

    fit.t.pred = fit.t.pred,
    cases = cases,
    scenario = scenario,
    H.E = H.E,
    H.W = H.W,
    adj.period = adj.period,
    population = population,
    rho = rho,
    serial_mean = serial_mean,
    serial_var = serial_var,
    lt = lt,
    window_size = window_size
  )

  yt.pr[ranges] <- abs(yt.pred$`Predicted Infections`)[ranges]
  rt.pr[ranges] <- yt.pred$`Predicted Rt`[ranges]
  lambdat.pr[ranges] <- yt.pred$`Predicted Lambda t`[ranges]

  wave.starts <-
    c(1, which(diff(ranges) > 1) + 1) # index of wave starts in ranges
  wave.ends <-
    c(which(diff(ranges) > 1), length(ranges))
  ## yt <- yt*rho
  ## yt.pr <- yt.pr*rho
  fit <- 0
  fit1 <- 0
  for (i in 1:length(wave.starts)) {
    curr.w.r <-
      ranges[wave.starts[i]:wave.ends[i]] # range of the current wave

    # Uses normal distribution for yt

    fit <- fit - sum(dnorm(
      x = yt[curr.w.r],
      mean = yt.pr[curr.w.r],
      sd = sqrt(yt.pr[curr.w.r] * abs(param[14])),
      log = TRUE
    ), na.rm = silence.errors)
  }

  return(fit)
}

#' @title Fit the Model
#' @description
#' Estimate the parameters of the model by maximizing the likelihood function (or, rather, by minimizing the negative log likelihood).
#'
#' @examples
#' library(REffectivePred)
#' ## Read in the data
#' path_to_data <- system.file("extdata/NY_OCT_4_2022.csv", package = "REffectivePred")
#' data <- read.csv(path_to_data)
#' head(data)
#' cases <- diff(c(0, data$cases)) # Convert cumulative cases into daily cases
#' lt <- length(cases)             # Length of cases
#' Time <- as.Date(data$date, tryFormats = c("%d-%m-%Y", "%d/%m/%Y"))
#'
#' navigate_to_config() # Open the config file, make any necessary changes here.
#' path_to_config <- system.file("config.yml", package = "REffectivePred")  # Read config file
#' cfg <- load_config()    # Build the cfg object
#'
#' ##### Option 1: populate the global environment with args to pass to function.
#' population <- cfg$population # Population size
#' window_size <- cfg$window.size
#' adj.period <- cfg$adj.period
#' fit.t.pred <- cfg$fit.t.pred # Time of prediction
#' not.predict <- cfg$not.predict
#' rt.func.num <- cfg$rt.func.num # choose which Rt function you want to use
#' num.iter <- cfg$num.iter
#' silence.errors <- cfg$silence.errors
#' predict.beyond <- cfg$predict.beyond
#' curve_params <- as.double(unlist(cfg$curve_params))
#' vt_params <- as.double(unlist(cfg$vt_params)) # The vt initial values, starting at wave 2
#' restriction_levels <- as.double(unlist(cfg$restriction_levels)) # Psi, u, and v parameters
#' betas <- as.double(unlist(cfg$betas)) #   betas
#' ini_params <- c(curve_params, vt_params, restriction_levels, betas)
#' restrictions_params <- cfg$restrictions_params
#' restriction_st_params <- cfg$restriction_st_params
#' param_scale <- abs(ini_params) / 10
#' waves_list <- ranges_to_waves(cfg$waves_list)
#' params_limits <- cfg$params_limits
#' num_waves <- cfg$num_waves
#' waves <- waves_1d_list(num_waves, waves_list)
#' rho <- eval(parse(text = cfg$rho))
#' serial_mean <- cfg$serial_mean
#' serial_var <- cfg$serial_var
#'
#' est <- estimate.mle(
#'   ini_params = ini_params,
#'   params_limits = params_limits,
#'   restrictions = restrictions_params,
#'   restriction.starts = restriction_st_params,
#'   ranges = waves,
#'   rt_func = rt.func.num,
#'   silence.errors = silence.errors,
#'
#'   fit.t.pred = fit.t.pred,
#'   param_scale = param_scale,
#'   num.iter = num.iter,
#'   cases = cases,
#'   scenario = NULL,
#'   H.E = NULL,
#'   H.W = NULL,
#'   adj.period = adj.period,
#'   population = population,
#'   rho = rho,
#'   serial_mean = serial_mean,
#'   serial_var = serial_var,
#'   lt = lt,
#'   window_size = window_size,
#'   hessian = FALSE
#' )
#' print(est)
#'
#' ##### Option 2: pass the cfg object instead.
#' est <- estimate.mle(
#'     cases = cases,
#'     cfg = cfg,
#'     hessian = FALSE
#'     )
#' print(est)
#' @param hessian Logical. If TRUE, computes the variance-covariance matrix at the MLE.
#' @param H.E Mobility metrics for category Retail & Entertainment. Currently unsupported.
#' @param H.W Mobility metrics for category Workplaces. Currently unsupported.
#' @param cases Vector of case counts.
#' @param cfg The object that contains all variables from the configuration file.
#'   This includes all function arguments except for \code{cases}, \code{hessian}, \code{H.E}, and \code{H.W}. All other
#'   arguments are overridden if \code{cfg} is passed to the method.
#' @param ini_params  Initial parameter values to be used in optimization. Includes the following sets of parameters in a vector, in this order:
#'  \itemize{
#'    \item   (a1,a2,a3,a4) = parameters for curve c() specified by \code{rt_func};
#'    \item   nu = loss of immunity rate;
#'    \item   (v2,v3,v4,v5) = transmissibility of variants in waves 2+, as relative multiplication factors compared to transmissibility in wave 1;
#'    \item   (psi1,psi2,psi3,psi4) = psi parameters for severity levels 1,2,3 and 4.
#'    \item   (u,v) = variance parameters. Only u is currently in use.
#'    \item   (beta0,beta.R,beta.E,beta.W), when restrictions = NULL. Currently unsupported.
#'    }
#' @param params_limits Boundaries/limits of the ini_params.
#' @param restrictions A numeric integer vector giving the severity of restrictions.
#'   Zero means no restriction, and higher numbers means greater severity/disruption.
#'   The ordered unique values should be consecutive integers starting from zero.
#'   Each number (other than 0) adds a new parameter to the fit. restrictions = NULL
#'   causes the function to use mobility data instead of the psi values (currently unsupported).
#' @param restriction.starts A vector of same length as restrictions, of times when restrictions
#'   came into effect. Note: the first index time should be 1.
#' @param ranges An vector of time ranges for the different waves.
#'   The waves ranges should be contiguous, with at least one unit of time
#'   between consecutive waves.
#' @param rt_func The parametric form of function c(). Options are listed under function c_helper.
#' @param silence.errors Logical. If TRUE, ignores certain errors to allow optimization to proceed. Not all errors can be ignored.
#' @param fit.t.pred Time from which prediction is done. If use.actual.not.predicted is TRUE, values of \eqn{S_t} before this time will be computed using actual counts.
#' @param param_scale Parameter scale. Passed as argument parscale to optim.
#' @param num.iter Maximum number of iterations. Passed as argument maxit to optim.
#' @param scenario A character string describing options to deal with restrictions. Currently unsupported.
#' @param adj.period Delays in society adjusting.
#' @param population total population size.
#' @param rho A vector of under-reporting rates of the same length as cases. If a scalar is supplied, the vector will be constant with this value.
#' @param serial_mean Mean of the serial interval on the log scale.
#' @param serial_var Variance of the serial interval on the log scale.
#' @param lt The length of cases.
#' @param window_size The maximum value for the serial interval.
#' @param verbose Logical. If TRUE, provides additional details while running the function.
#' @returns A list of maximum likelihood estimates of the parameters. Includes:
#'  \itemize{
#'    \item a1
#'    \item a2
#'    \item a3
#'    \item a4
#'    \item nu
#'    \item vt_params_est
#'    \item Psi
#'    \item betas
#'    \item negative_log_lik
#'    \item mle
#'    \item hessian
#'    \item SE
#'  }
#' @export
estimate.mle <- function(hessian = FALSE,
                         H.E = NULL,
                         H.W = NULL,
                         cases = NULL,

                         # Config object
                         cfg = NULL,

                         # None config object
                         ini_params = NULL,
                         params_limits = NULL,
                         restrictions = NULL,
                         restriction.starts = NULL,
                         ranges = NULL,
                         rt_func = 1,
                         silence.errors = FALSE,
                         fit.t.pred = NULL,
                         param_scale = NULL,
                         num.iter = NULL,
                         scenario = NULL,
                         adj.period = NULL,
                         population = NULL,
                         rho = NULL,
                         serial_mean = serial_mean,
                         serial_var = serial_var,
                         lt = NULL,
                         window_size = NULL,
                         verbose = FALSE) {

  if (!is.null(cfg)){
    population <- cfg$population # Population size

    lt <- length(cases) # Length of cases
    window_size <- cfg$window.size # Window size for Instantaneous Rt
    adj.period <- cfg$adj.period # Assume that psi shifts linearly to the new value over this period, to account for delays in society adjusting, and reporting delays

    fit.t.pred <- cfg$fit.t.pred #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Time of prediction
    rt_func <- cfg$rt.func.num # choose which Rt function you want to use
    num.iter <- cfg$num.iter
    silence.errors <- cfg$silence.errors
    predict.beyond <- cfg$predict.beyond

    curve_params <- as.double(unlist(cfg$curve_params))

    vt_params <- as.double(unlist(cfg$vt_params)) # The vt initial values, starting at wave 2
    restriction_levels <- as.double(unlist(cfg$restriction_levels)) #   Psi, u, and v parameters

    ini_betas <- as.double(unlist(cfg$betas)) #   betas

    ini_params <- c(curve_params, vt_params, restriction_levels, ini_betas)

    restrictions <- cfg$restrictions_params
    restriction.starts <- cfg$restriction_st_params
    param_scale <- abs(ini_params) / 10
    waves_list <- ranges_to_waves(cfg$waves_list)
    params_limits <- cfg$params_limits
    num_waves <- cfg$num_waves
    waves <- waves_1d_list(num_waves, waves_list)
    ranges <- waves
    rho <- eval(parse(text = cfg$rho)) # Under-reporting fraction (ref Irons, Raftery, 2021)
    serial_mean <- cfg$serial_mean
    serial_var <- cfg$serial_var

    infections <- cases / rho # Estimated time series of total infections
  }

  if (max(ranges) > fit.t.pred) {
    try(if(which(ranges > fit.t.pred)[1] < 1) stop("Error. fit.t.pred is before ranges starts."))
    ranges <- ranges[1:(which(ranges > fit.t.pred)[1] - 1)]
  }

  op <- optim(
    par = ini_params,
    fn = log_lklh,
    gr = NULL,
    fit.t.pred = fit.t.pred,
    rt_func = rt_func,
    params_limits = params_limits,
    restrictions = restrictions,
    restriction.starts = restriction.starts,
    ranges = ranges,
    silence.errors = silence.errors,

    lt = lt,
    cases = cases,
    scenario = scenario,
    H.E,
    H.W,
    adj.period = adj.period,
    population = population,
    rho = rho,
    serial_mean = serial_mean,
    serial_var = serial_var,
    window_size = window_size,
    method = "BFGS",
    control = list(
      trace = 6,
      parscale = param_scale,
      maxit = num.iter
    ),
    hessian = hessian
  )

  if(verbose){
    print(cat("Negative loglik of fit =", op$value, "\n"))
  }

  params <- modif.helper(op$par, params_limits)
  a1 <- params[1]
  a2 <- params[2]
  a3 <- params[3]
  a4 <- params[4]
  nu <- params[5]
  vt_params_est <- params[6:9]
  psi <- params[10:15]
  betas <- params[16:19]

  output <- list(
    "a1" = a1,
    "a2" = a2,
    "a3" = a3,
    "a4" = a4,
    "nu" = nu,
    "vt_params_est" = vt_params_est,
    "Psi" = psi,
    "betas" = betas,
    "negative_log_lik" = op$value,
    "mle" = params
  )

  if (hessian) {
    H <- op$hessian
    # Note: this is the hessian of the full likelihood, which will be like n*l(theta)
    #       where l is for a single obs; however, the final variance gets divided by n
    #       for iid data since sqrt(n)*(theta^ - theta) ~ N(0,I(theta)^-1). So nevermind the n.
    missing.params <- which(rowSums(H * H) == 0)
    missing.params <- which(rowSums(H * H) < 1e-10) # Added to solve a problem with values close to zero for nu for WPG.
    H.red <- H[-missing.params, -missing.params]
    varcov.red <- solve(H.red)
    se.vec <- rep(NA, length(op$par))
    se.vec[-missing.params] <- sqrt(diag(varcov.red))
    output <- c(output, list("hessian" = H, "SE" = se.vec))
  }
  return(output)
}

################## Prediction from real data##############
#' @title Epidemic Curve Model
#' @description
#' Computes the epidemic curve and associated quantities for a given parameter set.
#'
#' @examples
#' library(REffectivePred)
#' ## Read in the data
#' path_to_data <- system.file("extdata/NY_OCT_4_2022.csv", package = "REffectivePred")
#' data <- read.csv(path_to_data)
#' head(data)
#' cases <- diff(c(0, data$cases)) # Convert cumulative cases into daily cases
#' lt <- length(cases)             # Length of cases
#' Time <- as.Date(data$date, tryFormats = c("%d-%m-%Y", "%d/%m/%Y"))
#'
#' navigate_to_config() # Open the config file, make any necessary changes here.
#' path_to_config <- system.file("config.yml", package = "REffectivePred")  # Read config file
#' cfg <- load_config()    # Build the cfg object
#'
#' # Example 1. Using fits from Romanescu et al. (2023)
#'
#' r1 <- pred.curve(
#' a1 = 0.58,
#' a2 = 1.12,
#' nu = 0.56,
#' variant.transm = c(1,1.22,0.36,0.56),
#' Psi = c(0.58,0.52,0.49),
#' cases = cases,
#' cfg = cfg
#' )
#'
#' plot(cases, xlab="Day", ylab="Predicted cases")
#' lines(r1$'Predicted Cases', col='red')
#'
#' # Example 2. Best fit curve
#' est <- estimate.mle(
#'     cases = cases,
#'     cfg = cfg
#'     )
#' a1 <- est$a1
#' a2 <- est$a2
#' a3 <- est$a3
#' a4 <- est$a4
#' nu <- est$nu
#' vt <- c(1, est$vt_params_est)
#' psi <- est$Psi
#' betas <- est$betas
#'
#' r1 <- pred.curve(
#' a1 = a1,
#' a2 = a2,
#' a3 = a3,
#' a4 = a4,
#' nu = nu,
#' variant.transm = vt,
#' Psi = psi,
#' betas = betas,
#' cases = cases,
#' cfg = cfg
#' )
#' plot(r1$'Predicted Infections', xlab="Day", ylab="Predicted infections")
#'
#' @param a1,a2,a3,a4 Parameters of the contact rate curve specified by \code{rt_func}.
#' @param nu Loss of immunity rate beyond the first wave.
#' @param variant.transm Vector of transmissibility of variants in each wave,
#'   as relative multiplication factors compared to transmissibility in wave 1. Should always be 1 for the first wave.
#' @param Psi Vector of restriction parameters for severity levels 1 - 4.
#' @param betas Vector containing (beta0,beta.R,beta.E,beta.W), when restrictions = NULL. Not currently implemented.
#' @param cases vector of case counts.
#' @param cfg The object that contains all variables from the configuration file.
#'   \code{a1}, \code{a2}, \code{a3}, \code{a4}, \code{nu}, \code{variant.transm},
#'   \code{Psi}, \code{betas}, and \code{cases} are also required for the method
#'   to execute. All other parameters will not be used if \code{cfg} is passed
#'   to the method.
#' @param use.actual.not.predicted Logical; if FALSE (default), the susceptible fraction is updated using predicted cases. Otherwise updated using actual cases.
#' @param restrictions A numeric integer vector giving the severity of restrictions.
#'   Zero means no restriction, and higher numbers means greater severity/disruption.
#'   The ordered unique values should be consecutive integers starting from zero.
#'   Each number (other than 0) adds a new parameter to the fit.
#' @param restriction.starts A vector of same length as restrictions, of times when restrictions
#'   came into effect. Note: the first index time should be 1.
#' @param ranges A vector of time ranges for the different waves.
#'   The wave ranges should be contiguous, with at least one unit of time
#'   between consecutive waves.
#' @param rt_func The parametric form of function c(). Options are listed under function \code{c_helper}.
#' @param fit.t.pred Time from which prediction is done. If use.actual.not.predicted is TRUE, values of \eqn{S_t} before this time will be computed using actual counts.
#' @param predict.beyond Number of days to predict beyond the end of \code{cases}. See Details for usage notes.
#' @param scenario A character string describing options to deal with restrictions. Currently unsupported.
#' @param H.E Mobility metrics for category Retail & Entertainment. Currently unsupported.
#' @param H.W Mobility metrics for category Workplaces. Currently unsupported.
#' @param adj.period Adjustment period following a change in severity level. Restriction level (psi) is linearly interpolated from the old to the new value over this period.
#' @param population Total population size.
#' @param rho A vector of under-reporting rates of the same length as cases. If a scalar is supplied, the vector will be constant with this value.
#' @param serial_mean Mean of the serial interval on the log scale.
#' @param serial_var Variance of the serial interval on the log scale.
#' @param lt The length of cases.
#' @param window_size The maximum value for the serial interval.
#' @param verbose Logical. If TRUE, provides additional details while running the function.
#' @details At each time step, \eqn{R_{t}} is computed using the contact rate function \eqn{c(S_t)} implemented via \link{c_helper}. Then the number of cases is estimated using formula:
#'     \deqn{y_{t+1}=R_{t} \sum_{s=1}^M w_s y_{t+1-s}}
#'     Finally, the fraction \eqn{S_{t+1}} is updated. This creates a curve over the entire range of \code{ranges}. See Romanescu R, Hu S, Nanton D, Torabi M, Tremblay-Savard O, Haque MA. The effective reproductive number:
#'     modeling and prediction with application to the multi-wave Covid-19 pandemic. Epidemics. 2023 Jul 20:100708 \doi{10.1016/j.epidem.2023.100708} for more details.
#'
#'     For predicting an ongoing wave beyond the end of cases, the end of \code{ranges} (or \code{waves_list}, if using \code{cfg})
#'     should be specified to match the \code{predict.beyond} argument. As well, argument \code{use.actual.not.predicted} should be set to FALSE when predicting beyond the end of \code{cases}.
#' @return Returns list:
#'    \itemize{
#'      \item Predicted Infections - Vector of estimated infections, computed as predicted cases divided by rho.
#'      \item Predicted Cases - Vector of predicted cases.
#'      \item Predicted \eqn{R_t} - Vector of predicted susceptible fractions
#'      \item Predicted \eqn{R_t} - Vector of (model) predicted \eqn{R_t}.
#'      \item Predicted Lambda t - Vector of predicted Lambda_t, which is the numerator used in computing the empirical \eqn{R_t}.
#'      \item Psi.vec - Vector of psi values, which pastes together parameters psi over the period they apply, or 1 when there are no restrictions.
#'      \item Contact rate params - Vector of the curve parameters (a1, a2, a3, a4).
#'      }
#' @export
pred.curve <- function(a1 = 0,
                       a2 = 0,
                       a3 = 0,
                       a4 = 0,
                       nu = 0,
                       variant.transm = NULL,
                       Psi = NULL,
                       betas = NULL,
                       cases = NULL,

                       # With config object
                       cfg = NULL,

                       # No config object
                       use.actual.not.predicted = FALSE,
                       restrictions = NULL,
                       restriction.starts = NULL,
                       ranges = NULL,
                       rt_func = 1,
                       fit.t.pred = NULL,
                       predict.beyond = 0,
                       scenario = NULL,
                       H.E = NULL,
                       H.W = NULL,
                       adj.period = NULL,
                       population = NULL,
                       rho = NULL,
                       serial_mean = NULL,
                       serial_var = NULL,
                       lt = NULL,
                       window_size = NULL,
                       verbose = FALSE) {

  if (!is.null(cfg)){
    population <- cfg$population # Population size
    predict.beyond <- cfg$predict.beyond
    window_size <- cfg$window.size
    adj.period <- cfg$adj.period # Assume that psi shifts linearly to the new value over this period, to account for delays in society adjusting, and reporting delays

    fit.t.pred <- cfg$fit.t.pred #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Time of prediction
    use.actual.not.predicted <- cfg$not.predict
    rt_func <- cfg$rt.func.num # choose which Rt function to use
    num.iter <- cfg$num.iter
    silence.errors <- cfg$silence.errors
    predict.beyond <- cfg$predict.beyond

    curve_params <- as.double(unlist(cfg$curve_params))

    vt_params <- as.double(unlist(cfg$vt_params)) # The vt initial values, starting at wave 2
    # restriction_levels <- as.double(unlist(cfg$restriction_levels)) #   Psi (vector), u, v

    betas <- as.double(unlist(cfg$betas)) #   betas

    restrictions <- cfg$restrictions_params
    restriction.starts <- cfg$restriction_st_params

    params_limits <- cfg$params_limits
    num_waves <- cfg$num_waves # This doesn't get used in any function.
    waves_list <- ranges_to_waves(cfg$waves_list)
    waves <- waves_1d_list(num_waves, waves_list)
    ranges <- waves
    rho <- eval(parse(text = cfg$rho)) # Under-reporting fraction
    serial_mean <- cfg$serial_mean
    serial_var <- cfg$serial_var
    infections <- cases / rho # Estimated time series of total infections

    if (max(ranges) > fit.t.pred && FALSE) {
      try(if(which(ranges > fit.t.pred)[1] < 1) stop("Error. fit.t.pred is before ranges starts."))
      ranges <- ranges[1:(which(ranges > fit.t.pred)[1] - 1)]
    }
  }

  if (predict.beyond > 0) { # Note: this is needed if trying to predict past the last date.
    cases <- c(cases, rep(1, predict.beyond)) # Adding days to be able to predict
    if (length(cases) != ranges[length(ranges)]){
      warning("Warning: predict.beyond should match the end of the last wave")
    }
  }

  lt <- length(cases) # Length of cases

  if (length(rho) > 1 &&  length(rho) < lt) {
    rho <- c(rho, rep(rho[length(rho)], lt - length(rho)))  # Pad with last value if too short
  }

  if (length(rho) == 1) {
    rho <- rep(rho, lt)
  }

  infections <- cases / rho   # Estimated time series of total infections

  wave.starts <-
    c(1, which(diff(ranges) > 1) + 1) # *index* of wave starts in ranges
  wave.ends <- c(which(diff(ranges) > 1), length(ranges))

  # Define a vector of restriction severities, as long as yt.
  try(if(restriction.starts[1] != 1) stop("Error: restriction.starts should start at time 1."))

  yt.pr <- rep(NA, lt)
  rt.pr <- NA * yt.pr
  # This will hold the infections between waves, and the predicted yt within waves.
  st.pr <- NA * yt.pr
  ct.pr <- NA * yt.pr
  lambda.t.pr <- NA * yt.pr
  severities <- NA * yt.pr
  psi.vec.pub.health <- rep(1, length(yt.pr)) # holds the psi parameters (non-GMR related)
  psi.vec <- rep(1, length(yt.pr)) # holds the final/combined psi_t parameters for each index
  delta.psi.vec <- rep(0, length(yt.pr)) # psi[t] = psi[t-1] + delta.psi[t]


  for (i in 1:length(restrictions)) {
    severities[restriction.starts[i]:length(yt.pr)] <- restrictions[i]

    if (restrictions[i] > 0) {
      psi.vec.pub.health[restriction.starts[i]:length(yt.pr)] <- Psi[restrictions[i]]
      if (verbose){
        sprintf("%d", length(yt.pr))
      }
    } else {
      psi.vec.pub.health[restriction.starts[i]:length(yt.pr)] <- 1
    }

    # adjustment period
    if (i > 1) {
      st <- restriction.starts[i]
      en <- restriction.starts[i] + adj.period - 1
      psi.vec.pub.health[st:en] <-
        psi.vec.pub.health[st - 1] + (psi.vec.pub.health[st] - psi.vec.pub.health[st - 1]) / adj.period * (1:adj.period)
    }
  }

  psi.vec <- psi.vec.pub.health

  if (is.na(betas[1])) { # when no mobility data is available
    # will just use public health restrictions params
  } else { # Use Google data
    beta0 <- betas[1]
    beta.R <- betas[2]
    beta.E <- betas[3]
    beta.W <- betas[4]
  }

  for (i in 1:length(wave.starts)) { # For each wave

    if (i != length(wave.starts)){
      curr.w.r <-
        ranges[wave.starts[i]:wave.ends[i]] # range of the current wave
    }else{
      curr.w.r <-
        ranges[wave.starts[i]:(wave.ends[i] - 1)] # range of the current wave
    }

    # Prime (initialize) st.pr before waves & in-between waves
    if (i == 1) {
      # First wave
      st.pr[1:curr.w.r[1]] <-
        1 - cumsum(infections[1:curr.w.r[1]] / (population * 1)) # Before: * rho[tw + 1]
      yt.pr[1:curr.w.r[1]] <- infections[1:curr.w.r[1]]
      ct.pr[curr.w.r[1]] <- infections[curr.w.r[1]] / (population * 1) # Before: * rho[tw + 1]
      # For SIRS, Ct is reset to zero at the start of each wave, b/c the carrying capacity is St at said time.
    } else {
      # Subsequent waves
      st.pr[(ranges[wave.ends[i - 1]] + 1):(curr.w.r[1])] <-
        st.pr[ranges[wave.ends[i - 1]]] - cumsum(infections[(ranges[wave.ends[i - 1]] + 1):(curr.w.r[1])] / (population * 1)) # Before: * rho[tw + 1]
      yt.pr[(ranges[wave.ends[i - 1]] + 1):(curr.w.r[1])] <-
        infections[(ranges[wave.ends[i - 1]] + 1):(curr.w.r[1])]
      ct.pr[curr.w.r[1]] <- infections[curr.w.r[1]] / (population * 1) # Before: * rho[tw + 1]
    }


    if (i > 1 && nu > 0) {
      ## For waves 2+ with reinfection, loss of immunity is nu*(everyone immune or infected in the prev wave).
      rtA.pr <- rep(NA, lt)
      rtB.pr <- rtA.pr

      st.pr[(ranges[wave.ends[i - 1]] + 1):curr.w.r[1]] <-
        st.pr[ranges[wave.ends[i - 1]]] - cumsum(infections[(ranges[wave.ends[i - 1]] + 1):curr.w.r[1]] / (population * 1)) # Before: * rho[tw + 1]
      ST_im1 <- st.pr[curr.w.r[1] - 1]
      # Note: this computes st.pr between waves, assuming the decay dynamic continues until the beginning of the next wave.

      yt.pr[(ranges[wave.ends[i - 1]] + 1):curr.w.r[1]] <- infections[(ranges[wave.ends[i - 1]] + 1):curr.w.r[1]]

      # Updating st.pr just before a new wave with ppl who lost immunity
      st.pr[curr.w.r[1]] <-
        (1 - nu) * st.pr[curr.w.r[1]] + nu

      wA <- max(1 - nu / st.pr[curr.w.r[1]], 0)
    }

    for (tw in curr.w.r) {
      # At the very beginning after time 0.
      if (tw < (window_size)) {
        # for n=2:9 n=1 not possible due to unavailability of rt
        lambda.t.pr[tw] <- sum(c(rep(0, window_size - tw), yt.pr[1:tw]) * serial.helper(window_size = window_size, serial_mean = serial_mean, serial_var = serial_var))
      } else {
        lambda.t.pr[tw] <- sum(yt.pr[(tw - window_size + 1):tw] * serial.helper(window_size = window_size, serial_mean = serial_mean, serial_var = serial_var))
      }


      # ***************************** Module uses Google Mobility Reports


      if (!is.na(betas[1])) {
        #     solve for psi_t
        # This is a stub
      } # ***************************** End module


      # Update step without SIRS mixing assumptions. For waves 2+ it will be replaced below, using the correct ass'ns.
      rt.pr[tw] <- variant.transm[i] * st.pr[tw] * c_helper(
        rt_func = rt_func,
        st.inner = st.pr[tw],
        a1 = a1, a2 = a2, a3 = a3, a4 = a4,
        psi = psi.vec[tw]
      )
      rt.pr[rt.pr < 0] <- 0

      yt.pr[tw + 1] <-
        rt.pr[tw] * lambda.t.pr[tw]
      if (tw < fit.t.pred &&
        use.actual.not.predicted) {
        # replace with actual
        yt.pr[tw + 1] <- infections[tw + 1]
      }
      st.pr[tw + 1] <- st.pr[tw] - yt.pr[tw + 1] / (population * 1) # Before: * rho[tw + 1]
      ct.pr[tw + 1] <- ct.pr[tw] + yt.pr[tw + 1] / (population * 1) # Before: * rho[tw + 1]

      if (i > 1 && nu > 0 && !rt_func %in% c(3)) {
        # First day of new wave
        # This is taken care of before, in lambda.t.pr[tw] = ...

        rtA.pr[tw] <- variant.transm[i] * st.pr[tw] * c_helper(
          rt_func = rt_func,
          st.inner = (ST_im1 * st.pr[tw] / st.pr[curr.w.r[1]]),
          a1 = a1, a2 = a2, a3 = a3, a4 = a4,
          psi = psi.vec[tw]
        )
        rtB.pr[tw] <- variant.transm[i] * st.pr[tw] * c_helper(
          rt_func = rt_func,
          st.inner = st.pr[tw] / st.pr[curr.w.r[1]],
          a1 = a1, a2 = a2, a3 = a3, a4 = a4,
          psi = psi.vec[tw]
        )

        rtA.pr[rtA.pr < 0] <- 0
        rtB.pr[rtB.pr < 0] <- 0

        yt.pr[tw + 1] <-
          (wA * rtA.pr[tw] + (1 - wA) * rtB.pr[tw]) * lambda.t.pr[tw]

        if (tw < fit.t.pred &&
          use.actual.not.predicted) {
          # replace with actual
          yt.pr[tw + 1] <- infections[tw + 1]
        }

        st.pr[tw + 1] <- st.pr[tw] - yt.pr[tw + 1] / (population * 1)
      }
    } # End loop on tw

    # Compute a harmonized Rt.pr
    missing.rt <- which(is.na(rt.pr))
    missing.rt <- missing.rt[missing.rt < lt]
    rt.pr[missing.rt] <-
      yt.pr[missing.rt + 1] / lambda.t.pr[missing.rt]
  }

  output <- list(
    "Predicted Infections" = yt.pr[1:lt],
    # Yt is now estimated infections, not just observed cases.
    "Predicted Cases" = rho * yt.pr[1:lt],
    "Predicted St" = st.pr[1:length(st.pr)],
    # "Predicted Ct" = ct.pr[1:length(ct.pr)],
    "Predicted Rt" = rt.pr,
    "Predicted Lambda t" = lambda.t.pr,
    "Psi.vec.pub.health" = psi.vec.pub.health,
    "Psi.vec" = psi.vec,
    "Contact rate params" = c(a1, a2, a3, a4)
  )
  output
}

## Function to compute partial derivatives of pred.curve at each time point. Internal. Do not put in package manual.
# Note: will not capture any change related to sigma.
D.pred.yt.1 <- function(param,
                        restrictions = NULL,
                        restriction.starts = NULL,
                        ranges = NULL,
                        rt_func = 1,
                        fit.t.pred,
                        predict.beyond = 0,
                        lt,
                        cases,
                        scenario = NULL,
                        H.E,
                        H.W,
                        adj.period,
                        population,
                        rho,
                        serial_mean = serial_mean,
                        serial_var = serial_var,
                        window_size,
                        eps = .Machine$double.eps ^ (1/2),
                        deriv.order) {

  param0 <- param

  yt.pred <- pred.curve(
    a1 = param0[1],
    a2 = param0[2],
    a3 = param0[3],
    a4 = param0[4],
    nu = param0[5],
    Psi = param0[10:15],
    betas = param0[16:19],
    restrictions = restrictions,
    restriction.starts = restriction.starts,
    ranges = ranges,
    variant.transm = c(1, param0[6:9]) ,
    rt_func = rt_func,

    fit.t.pred = fit.t.pred,
    predict.beyond = predict.beyond,
    cases = cases,
    scenario = scenario,
    H.E = H.E,
    H.W = H.W,
    adj.period = adj.period,
    population = population,
    rho = rho,
    serial_mean = serial_mean,
    serial_var = serial_var,
    lt = lt,
    window_size = window_size
  )$"Predicted Cases"

  if (deriv.order == 0){
    return(yt.pred)
  }
  if (deriv.order == 1){
    D = matrix(NA, nrow=length(param), ncol=lt + predict.beyond)
    for (i in 1:length(param)){
      param.plus = param
      param.plus[i] = param[i] + eps

      d = (D.pred.yt.1(param.plus,
                       restrictions = restrictions,
                       restriction.starts = restriction.starts,
                       ranges = ranges,
                       rt_func = rt_func,
                       fit.t.pred = fit.t.pred,
                       predict.beyond = predict.beyond,
                       cases = cases,
                       scenario = scenario,
                       H.E = H.E,
                       H.W = H.W,
                       adj.period = adj.period,
                       population = population,
                       rho = rho,
                       serial_mean = serial_mean,
                       serial_var = serial_var,
                       lt = lt,
                       window_size = window_size,
                       deriv.order = 0) - yt.pred) / eps

      D[i,] = d
    }
    return(D)
  }
  try(if(!deriv.order %in% c(0,1)) stop("Unrecognized derivative order."))
}




#' @title Confidence bands
#' @description Computes the pointwise confidence interval of the epidemic curve.
#'
#' @examples
#' library(REffectivePred)
#' ## Read in the data
#' path_to_data <- system.file("extdata/NY_OCT_4_2022.csv", package = "REffectivePred")
#' data <- read.csv(path_to_data)
#' head(data)
#' cases <- diff(c(0, data$cases)) # Convert cumulative cases into daily cases
#' lt <- length(cases)             # Length of cases
#' Time <- as.Date(data$date, tryFormats = c("%d-%m-%Y", "%d/%m/%Y"))
#'
#' navigate_to_config() # Open the config file, make any necessary changes here.
#' path_to_config <- system.file("config.yml", package = "REffectivePred")  # Read config file
#' cfg <- load_config()    # Build the cfg object
#'
#' # Estimate parameters
#' est <- estimate.mle(
#'     cases = cases,
#'     cfg = cfg,
#'     hessian = TRUE
#'     )
#' a1 <- est$a1
#' a2 <- est$a2
#' a3 <- est$a3
#' a4 <- est$a4
#' nu <- est$nu
#' vt <- c(1, est$vt_params_est)
#' psi <- est$Psi
#' betas <- est$betas
#'
#' # Predict curve
#' r1 <- pred.curve(
#' a1 = a1,
#' a2 = a2,
#' a3 = a3,
#' a4 = a4,
#' nu = nu,
#' variant.transm = vt,
#' Psi = psi,
#' betas = betas,
#' cases = cases,
#' cfg = cfg
#' )
#'
#' plot_outputs(Time = Time,
#' cases = cases,
#' cfg = cfg,
#' curve = r1,
#' option = 2
#' )
#'
#' bounds <- ci.curve(fit = est,
#'                    cases = cases,
#'                    cfg = cfg)
#'
#' # Adding CI bands
#' # lines(c(Time, Time[length(Time)]+(1:predict.beyond)), bounds[2,], lty = 2)
#' # lines(c(Time, Time[length(Time)]+(1:predict.beyond)), bounds[1,], lty = 2)
#'
#'
#' @param fit Output from function estimate.mle.
#' @param H.E Mobility metrics for category Retail & Entertainment. Currently unsupported.
#' @param H.W Mobility metrics for category Workplaces. Currently unsupported.
#' @param scenario A character string describing options to deal with restrictions. Currently unsupported.
#' @param cases vector of case counts.
#' @param cfg The object that contains all variables from the configuration file.
#'   \code{fit}, \code{H.E}, \code{H.W}, \code{scenario}, and \code{cases} are also required
#'   for the method to execute. All other parameters will not be used if \code{cfg} is passed to the method.
#' @param restrictions A numeric integer vector giving the severity of restrictions.
#'   Zero means no restriction, and higher numbers means greater severity/disruption.
#'   The ordered unique values should be consecutive integers starting from zero.
#'   Each number (other than 0) adds a new parameter to the fit.
#' @param restriction.starts A vector of same length as restrictions, of times when restrictions
#'   came into effect. Note: the first index time should be 1.
#' @param ranges A vector of time ranges for the different waves.
#'   The wave ranges should be contiguous, with at least one unit of time
#'   between consecutive waves.
#' @param rt_func The parametric form of function c(). Options are listed under function c_helper.
#' @param predict.beyond Number of days to predict beyond the end of \code{cases}. See Details for usage notes.
#' @param fit.t.pred Time from which prediction is done. If use.actual.not.predicted is TRUE, values of \eqn{S_t} before this time will be computed using actual counts.
#' @param adj.period Adjustment period following a change in severity level. Restriction level (psi) is linearly interpolated from the old to the new value over this period.
#' @param population Total population size.
#' @param rho A vector of under-reporting rates of the same length as cases. If a scalar is supplied, the vector will be constant with this value.
#' @param serial_mean Mean of the serial interval on the log scale.
#' @param serial_var Variance of the serial interval on the log scale.
#' @param lt The length of cases.
#' @param window_size The maximum value for the serial interval.
#' @param eps The epsilon value for computing finite differences.
#' @return Returns a matrix with two rows containing Wald-style confidence bounds:
#'    \itemize{
#'      \item ci_lower - lower bound of a 95\% pointwise CI of the best fit curve.
#'      \item ci_upper - upper bound of a 95\% pointwise CI of the best fit curve.
#'    }
#' @export
ci.curve <- function(fit = NULL,
                     H.E = NULL,
                     H.W = NULL,
                     scenario = NULL,
                     cases = NULL,

                     cfg = NULL,

                     restrictions = NULL,
                     restriction.starts = NULL,
                     ranges = NULL,
                     rt_func = 1,
                     fit.t.pred = NULL,
                     predict.beyond = 0,
                     lt = NULL,
                     adj.period = NULL,
                     population = NULL,
                     rho = NULL,
                     serial_mean = NULL,
                     serial_var = NULL,
                     window_size = NULL,
                     eps = .Machine$double.eps ^ (1/2)){

  if (!is.null(cfg)){
    population <- cfg$population # Population size
    predict.beyond <- cfg$predict.beyond
    window_size <- cfg$window.size # Window size for Instantaneous Rt
    adj.period <- cfg$adj.period # Assume that psi shifts linearly to the new value over this period, to account for delays in society adjusting, and reporting delays

    fit.t.pred <- cfg$fit.t.pred #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Time of prediction
    rt_func <- cfg$rt.func.num # choose which Rt function you want to use

    restrictions <- cfg$restrictions_params
    restriction.starts <- cfg$restriction_st_params

    num_waves <- cfg$num_waves
    waves_list <- ranges_to_waves(cfg$waves_list)
    waves <- waves_1d_list(num_waves, waves_list)
    ranges <- waves
    rho <- eval(parse(text = cfg$rho))
    serial_mean <- cfg$serial_mean
    serial_var <- cfg$serial_var

    if (max(ranges) > fit.t.pred && FALSE) {  # Dec 17, 2023. Not sure that it is necessary to reduce ranges.
      try(if(which(ranges > fit.t.pred)[1] < 1) stop("Error. fit.t.pred is before ranges starts."))
      ranges <- ranges[1:(which(ranges > fit.t.pred)[1] - 1)]
    }
  }

  lt <- length(cases) # Length of cases

  params.mle <- fit$mle

  H <- fit$hessian
  D0 <- D.pred.yt.1(params.mle,
                    restrictions = restrictions,
                    restriction.starts = restriction.starts,
                    ranges = ranges,
                    rt_func = rt_func,
                    fit.t.pred = fit.t.pred,
                    predict.beyond = predict.beyond,
                    cases = cases,
                    scenario = scenario,
                    H.E = H.E,
                    H.W = H.W,
                    adj.period = adj.period,
                    population = population,
                    rho = rho,
                    serial_mean = serial_mean,
                    serial_var = serial_var,
                    lt = lt,
                    window_size = window_size,
                    eps = eps,
                    deriv.order = 0)
  D1 <- D.pred.yt.1(params.mle,
                    restrictions = restrictions,
                    restriction.starts = restriction.starts,
                    ranges = ranges,
                    rt_func = rt_func,
                    fit.t.pred = fit.t.pred,
                    predict.beyond = predict.beyond,
                    cases = cases,
                    scenario = scenario,
                    H.E = H.E,
                    H.W = H.W,
                    adj.period = adj.period,
                    population = population,
                    rho = rho,
                    serial_mean = serial_mean,
                    serial_var = serial_var,
                    lt = lt,
                    window_size = window_size,
                    eps = eps,
                    deriv.order = 1)
  D1 <- t(D1)  # this has #param columns, each of length lt

  missing.params = which(rowSums(H*H) == 0 )
  missing.params = which(rowSums(H*H) < 1e-10 )    # Added to solve a problem with values close to zero for nu.
  H.red = H[-missing.params, -missing.params]
  varcov.red <- solve(H.red)

  D1.red <- D1[,-missing.params]

  Var.ft <- D1.red %*% varcov.red %*% t(D1.red)  # For each time t, this conforms to dell(f_t) orientation in paper.

  sigmaP = sqrt(diag(Var.ft) + 0*fit$mle[14]*D0)
  ci_lower <- D0 - 1.96*sigmaP
  ci_upper <- D0 + 1.96*sigmaP
  ci <- data.frame(cbind(ci_lower, ci_upper))
  return(t(ci))
}


#' @title Contact rate function.
#' @description Computes the c() function.
#'
#' @param rt_func Options are:
#'  \itemize{
#'    \item 1 - Two piece exponential.
#'    \item 2 - Exponential power model adapted from Granich et al. (2009)
#'    \item 3 - Mass action.
#'    \item 4 - Shifted inverse.
#'    \item 5 - Power.
#'    \item 6 - Poisson.
#'    \item 7 - Geometric.
#'  }
#'
#' @param st.inner The susceptible fraction \eqn{S_t}.
#' @param a1,a2,a3,a4 Parameters of the contact rate curve specified by \code{rt_func};
#' @param psi A vector of same length as st.inner containing the corresponding psi restriction factor, or a scalar.
#' @details See Romanescu et al. (2023) for the exact forms of the functions.
#' @returns The value of the contact rate, used to compute \eqn{R_t}.
#' @export

c_helper <- function(rt_func = 1, st.inner = NULL,
                     a1 = NULL, a2 = NULL, a3 = NULL, a4 = NULL,
                     psi = NULL) {
  switch(rt_func,
    "1" = { # Two piece expl
      # Old version (2022) with 4 parameters: return((exp(a1 * (psi * st.inner - a2)) - 1 + max(exp(a3 * (psi * st.inner - a4)) - 1, 0)) )
      # New version (3 params)
      return(exp(a1 * (psi * st.inner - 0)) - 1 + max(exp(a3 * (psi * st.inner - a2)) - 1, 0))
    },
    "2" = { # Granich
      return(a1 * exp(-a2 * (1 - psi * st.inner)^a3))
    },
    "3" = {
      return(a1 * psi)
    },
    "4" = { # Shifted inverse
      return(a1 / pmax((a2 - psi * st.inner), 0.00001))
    },
    "5" = { # Power -- fit = 3,697 for NY, so better than MA but not great.
      return(a1 * (psi * st.inner)^a2)
    },
    "6" = { # Poisson (a2 = lambda)
      return(a1 * (a2 + log(psi * st.inner)))
    },
    "7" = { # Geometric
      return(2 * a1 * (psi * st.inner - 1 + exp(-1 / a2)) / (1 - exp(-1 / a2)))
    },
    stop("Invalid Rt function number.")
  )
}


######################### HELPER FUNCTION FOR modif ###########################
#' Helper function which ensures that parameters are within specified bounds. Called by the estimate.mle.
#' Note: this is an internal function which should not be modified.
#' @keywords internal
#' @param params_modif Raw parameter vector, in the same order as argument ini_param in estimate.mle.
#' @param params_limits Boundaries/limits of the ini_params.
#' @returns The modified parameters (bounded, if needed).

modif.helper <- function(params_modif, params_limits) {
  total_size <- length(params_limits)
  # restriction_levels_size <- length(restriction_levels)
  # rest_size <- total_size - restriction_levels_size
  # psi_positions = (rest_size + 1):(total_size - 2)
  psi_positions <- 10:13
  all_positions <- 1:total_size

  for (pos in all_positions[!all_positions %in% psi_positions]) {
    # absolute value is true
    if (params_limits[[pos]][[1]]) {
      params_modif[[pos]] <- abs(params_modif[[pos]])
    }

    # lower-bound exists
    if (is.numeric(params_limits[[pos]][[2]])) {
      params_modif[[pos]] <- max(params_limits[[pos]][[2]], params_modif[[pos]])
    }

    # upper-bound exists
    if (is.numeric(params_limits[[pos]][[3]])) {
      params_modif[[pos]] <- min(params_modif[[pos]], params_limits[[pos]][[3]])
    }
  }

  upperbound_tmp <- 1
  for (pos in psi_positions) {
    # absolute value is true
    if (params_limits[[pos]][[1]]) {
      params_modif[[pos]] <- abs(params_modif[[pos]])
    }

    # lower-bound exists
    if (is.numeric(params_limits[[pos]][[2]])) {
      # "upperbound_tmp" should be equal to the previous psi or 1 if it is the first loop.
      # When lowerbound of the current psi is larger than the upper_tmp, we might have a problem in limits.
      try(if(upperbound_tmp < params_limits[[pos]][[2]]) stop("Please reset boundaries."))

      params_modif[[pos]] <- max(params_limits[[pos]][[2]], params_modif[[pos]])
    }

    # upper-bound exists
    if (is.numeric(params_limits[[pos]][[3]])) {
      upperbound <- min(params_limits[[pos]][[3]], upperbound_tmp)
      params_modif[[pos]] <- min(params_modif[[pos]], upperbound)

      # set the upperbound_tmp
      upperbound_tmp <- params_modif[[pos]]
    }
  }

  return(params_modif)
}


#' @title Serial interval
#' @description Helper function for computing weights \code{w} based on the serial interval.
#'
#' @param window_size The maximum value for the serial interval.
#' @param serial_mean Mean of the serial interval on the log scale. See Details.
#' @param serial_var Variance of the serial interval on the log scale. See Details.
#' @details
#'  Computed based on a log normal density function, as in Nishiura et al. (2020). Parameters
#'  \code{serial_mean} and \code{serial_var} are arguments \code{meanlog} and \code{sdlog} of
#'  function \code{dlnorm}. Default values are taken from the same reference.
#'
#' @returns A vector that stores the serial interval in reverse order. This is
#'  meant to be multiplied to infections in chronological order.

serial.helper <- function(window_size, serial_mean = log(4.0), serial_var = log(1.380715)) {

  w <- dlnorm(window_size:1, meanlog = serial_mean, sdlog = sqrt(serial_var))
  return(w)
}

################################################################################


################################### PLOT FUNCTION#############################################


#' @title Plotting function
#' @description
#' Various plots related to an epidemic curve.
#'
#' @examples
#' library(REffectivePred)
#' ## Read in the data
#' path_to_data <- system.file("extdata/NY_OCT_4_2022.csv", package = "REffectivePred")
#' data <- read.csv(path_to_data)
#' head(data)
#' cases <- diff(c(0, data$cases)) # Convert cumulative cases into daily cases
#' lt <- length(cases)             # Length of cases
#' Time <- as.Date(data$date, tryFormats = c("%d-%m-%Y", "%d/%m/%Y"))
#'
#' navigate_to_config() # Open the config file, make any necessary changes here.
#' path_to_config <- system.file("config.yml", package = "REffectivePred")  # Read config file
#' cfg <- load_config()    # Build the cfg object
#'
#' # Estimate parameters
#' est <- estimate.mle(
#'     cases = cases,
#'     cfg = cfg
#'     )
#' a1 <- est$a1
#' a2 <- est$a2
#' a3 <- est$a3
#' a4 <- est$a4
#' nu <- est$nu
#' vt <- c(1, est$vt_params_est)
#' psi <- est$Psi
#' betas <- est$betas
#'
#' # Predict curve
#' r1 <- pred.curve(
#' a1 = a1,
#' a2 = a2,
#' a3 = a3,
#' a4 = a4,
#' nu = nu,
#' variant.transm = vt,
#' Psi = psi,
#' betas = betas,
#' cases = cases,
#' cfg = cfg
#' )
#'
#' plot_outputs(Time = Time,
#' cases = cases,
#' window_size = cfg$window.size,
#' serial_mean = cfg$serial_mean,
#' serial_var = cfg$serial_var,
#' predict.beyond = cfg$predict.beyond,
#' waves_list = cfg$waves_list,
#' num_waves = cfg$num_waves,
#' rt_func = cfg$rt.func.num,
#' curve = r1,
#' restrictions = cfg$restrictions_params,
#' restriction.starts = cfg$restriction_st_params,
#' rt.max = 10
#' )
#' @param curve The output list from the prediction function, see \code{pred.curve}.
#' @param Time A vector of dates corresponding to cases.
#' @param cases A vector containing cases for each time point.
#' @param cfg The object that contains all variables from the configuration file.
#'   \code{curve}, \code{Time}, and \code{cases} are also required
#'   for the method to execute. All other parameters will not be used if \code{cfg} is passed to the method.
#' @param window_size The maximum value for the serial interval.
#' @param serial_mean Mean of the serial interval on the log scale.
#' @param serial_var Variance of the serial interval on the log scale.
#' @param predict.beyond How many days to predict beyond the end of \code{cases}.
#' @param waves_list A two-dimensional list containing the waves' time data.
#' @param num_waves Total number of waves.
#' @param rt_func A flag that indicates which rt function to use. Should match the shape of \code{curve}.
#' @param restrictions A numeric integer vector giving the severity of restrictions.
#' @param restriction.starts A vector of same length as restrictions, of times when restrictions
#'   came into effect. Note: the first index time should be 1.
#' @param a1,a2,a3,a4 Parameters of the contact rate curve specified by \code{rt_func}. These override
#'  the values given in \code{curve} for the last plot only. If not specified, will use the values from \code{curve}.
#' @param rt.max An optional upper limit for the y-axis when plotting \eqn{R_t}.
#' @param option A choice of which plot to return (1,2, or 3 - see Value for options). If set to "all" (the default) plots all three figures.
#' @param verbose Logical. If TRUE, provides additional details while running the function.
#' @return NULL. Generates a few plots: a plot of \eqn{R_t} over time, with waves shaded (for \code{option} = 1); the epidemic curve
#' overlaid on top of observed cases (\code{option} = 2), where the shading reflects restriction measures; and a plot of the
#' theoretical \eqn{R_t} versus \eqn{S_t}, in a fully susceptible population with no restrictions (\code{option} = 3).
#' @export
plot_outputs <- function(curve = NULL,
                         Time = NULL,
                         cases = NULL,

                         cfg = NULL,

                         window_size = NULL,
                         serial_mean,
                         serial_var,
                         predict.beyond = 0,
                         waves_list = NULL,
                         num_waves = NULL,
                         rt_func = NULL,
                         restrictions = NULL,
                         restriction.starts = NULL,
                         a1 = NULL, a2 = NULL, a3 = NULL, a4 = NULL,
                         rt.max = NULL,
                         option = "all",
                         verbose = FALSE) {

  if (!is.null(cfg)){
    lt <- length(cases) # Length of cases
    window_size <- cfg$window.size # Window size for Instantaneous Rt
    serial_mean <- cfg$serial_mean
    serial_var <- cfg$serial_var
    rt_func <- cfg$rt.func.num # choose which Rt function you want to use

    restrictions <- cfg$restrictions_params
    restriction.starts <- cfg$restriction_st_params
    waves_list <- ranges_to_waves(cfg$waves_list)
    num_waves <- cfg$num_waves
    predict.beyond <- cfg$predict.beyond
  }

  # Estimate Rt using serial interval
  rt.empirical <- rt_empirical(cases, window_size, serial_mean, serial_var)
  # And predicted Rt
  Rt.pr <- curve$`Predicted Rt`

  lt <- length(curve$`Predicted Cases`)
  range.points <- 1:lt
  if (predict.beyond > 0){
    Time <- c(Time, Time[length(Time)] + (1:predict.beyond))
  }

  waves = waves_1d_list(num_waves, waves_list)

  # Plot predicted Rt (lines) on top of empirical Rt (points);
  if (is.null(rt.max)){   rt.max = max(rt.empirical, na.rm = TRUE)  }

  if (option == 1 || option == "all"){
    plot(Time[range.points], rt.empirical[range.points],
         xlab = "Time", ylab = "Rt", main = expression("Empirical and model R"[t]), # xlim=c(0, min(fit.t.pred+60,lt)),
         type = "n", ylim = c(0, rt.max), lwd = 2, xaxt = "n", cex.lab = 1.3, col = "darkgrey"
    )
    axis(1, Time[range.points], format(Time[range.points], "%b-%Y"), cex.axis = .9, tick = FALSE)

    # Shade in the waves.
    for (i in 1:num_waves) {
      rect(Time[min(waves_list[[i]])], 0, Time[max(waves_list[[i]])], rt.max, col = "mistyrose", lty = 0)
    }


    points(Time[range.points], rt.empirical[range.points], lwd = 2, col = "darkgrey", type = "p")

    lines(Time[range.points], Rt.pr[range.points], type = "l", xlab = "", ylab = "", lty = 1, col = "darkblue", lwd = 3, bty = "n") # First n points are observed
  }

  # Plot cases
  yt.pr1 <- abs(curve$`Predicted Cases`)
  yt.pr1.plot <- rep(NA, lt)
  yt.pr1.plot[waves] <- yt.pr1[waves]

  if (option == 2 || option == "all"){
    plot(Time[range.points], NA * cases[range.points],
         type = "n", ylim = c(0, max(cases) + 0 / 2),
         xlab = "Time", ylab = "Daily cases", main = "Epidemic curve", lwd = 2, xaxt = "n", cex.lab = 1.3
    )

    if (TRUE) {
      # Shade in the public health restrictions.
      colfunc <- colorRampPalette(c("azure1", "lightcyan3"))
      for (i in 1:(length(restrictions) - 1)) {
        if (restrictions[i] > 0) {
          rect(Time[restriction.starts[i]], 0, Time[restriction.starts[i + 1]], 1.3 * max(cases),
               col = colfunc(max(restrictions))[restrictions[i]], lty = 0
          )
        }
      }
    }
    points(Time[range.points], cases[range.points], xaxt = "n", col = "gray55", cex.lab = 1.3)
    axis(1, Time, format(Time, "%b-%Y"), cex.axis = .9, tick = FALSE)

    lines(Time[range.points], yt.pr1.plot[range.points], type = "l", xlab = "", ylab = "", lty = 1, col = "darkblue", lwd = 3, bty = "n") # First n points are observed
    if (verbose){
      print(yt.pr1.plot[range.points]) # test
    }
  }


  # Plot theoretical Rt in a fully susceptible population, with no restrictions
  ct.pr.plot <- (1:1000) / 1000
  st.pr.plot <- 1 - ct.pr.plot

  contact.rate.params <- curve$`Contact rate params`
  if (is.null(a1)){    a1 <- contact.rate.params[[1]] }
  if (is.null(a2)){    a2 <- contact.rate.params[[2]] }
  if (is.null(a3)){    a3 <- contact.rate.params[[3]] }
  if (is.null(a4)){    a4 <- contact.rate.params[[4]] }

  psi.scalar <- 1
  rt.pr.x.st <- st.pr.plot * c_helper(
    rt_func = rt_func,
    st.inner = st.pr.plot,
    a1 = a1, a2 = a2, a3 = a3, a4 = a4,
    psi = psi.scalar
  )

  if (option == 3 || option == "all"){
    rt.pr.plot <- rt.pr.x.st
    rt.pr.plot[rt.pr.plot < 0] <- NA
    plot(st.pr.plot, rt.pr.plot,
         type = "l", lwd = 2,
         col = "blue", xlab = "St", ylab = "Rt", xlim = c(0 + 0 * min(curve$`Predicted St`, na.rm = TRUE), 1), ylim = c(0, max(rt.pr.plot, na.rm = TRUE))
    )
    title(expression("Theoretical R"[t](S[t])))
    abline(h = 1, lty = 2, col = "purple")
  }
}
