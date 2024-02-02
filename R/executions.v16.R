##################################### MAIN #####################################
# includes underreporting estimation
#' @import yaml
#' @import config
#' @importFrom zoo as.Date
#' @importFrom grDevices dev.off
#' @importFrom utils data read.csv
#' @include prediction_model.v16.R
NULL

#' @title Navigate to the config file
#' @description  Prints the path to the config file and opens the config file.
#' @returns The path to the configuration file.
#' @export
navigate_to_config <- function() {
  path_to_config <- system.file("config.yml", package = "REffectivePred")
  system(path_to_config)
  output <- path_to_config
}

#' @title Load configuration file
#' @description Load the configuration file as an object in the global
#'    environment. This can be passed to the main functions instead of the
#'    individual arguments within, for user convenience.
#' @returns A ’config’ object, which is a list that stores input parameters and
#'    settings from config.yml. Variable names are imported exactly. See the
#'    configuration file for details.
#' @export
load_config <- function(){
  path_to_config <- system.file("config.yml", package = "REffectivePred")
  cfg <- config::get(file = path_to_config)
  output <- cfg
}


#' @title Demo of main functions
#' @description Fits the model to an example case data and predicts the epidemic curves and plots the outputs.
#' @details
#'  Please modify the config file before invoking this method.
#'  The config file contains all the settings and initial  parameter values
#'  necessary for the algorithm to run. Path to the dataset (in csv format)
#'  is also set in the config file. This file should be updated to the desired
#'  specifications before running this demo. To open it, execute
#'  REffectivePred::navigate_to_config().
#' @usage NULL
#' @param path_to_data Absolute path to the dataset in csv format.
#' @returns No return value.
#' @export
re_predict <- function(path_to_data = NULL) {
  path_to_config <- system.file("config.yml", package = "REffectivePred")
  cfg <- config::get(file = path_to_config)

  if (is.null(path_to_data)) {
    path_to_data <- system.file("extdata/NY_OCT_4_2022.csv", package = "REffectivePred")
  }

  data <- read.csv(path_to_data)
  cases <- diff(c(0, data$cases)) # Convert cumulative cases into daily cases
  Time <- as.Date(data$date, tryFormats = c("%d-%m-%Y", "%d/%m/%Y"))
  population <- cfg$population # Population size

  lt <- length(cases) # Length of cases
  window_size <- cfg$window.size # Window size for Instantaneous Rt
  rho <- eval(parse(text = cfg$rho)) # Under-reporting fraction (ref Irons, Raftery, 2021)
  serial_mean <- cfg$serial_mean
  serial_var <- cfg$serial_var
  adj.period <- cfg$adj.period # Assume that psi shifts linearly to the new value over this period, to account for delays in society adjusting, and reporting delays

  fit.t.pred <- cfg$fit.t.pred #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Time of prediction
  not.predict <- cfg$not.predict
  rt.func.num <- cfg$rt.func.num
  num.iter <- cfg$num.iter
  silence.errors <- cfg$silence.errors
  predict.beyond <- cfg$predict.beyond

  curve_params <- as.double(unlist(cfg$curve_params))

  vt_params <- as.double(unlist(cfg$vt_params)) # The vt initial values, starting at wave 2

  restriction_levels <- as.double(unlist(cfg$restriction_levels)) #   Psi

  betas <- as.double(unlist(cfg$betas)) #   betas

  ini_params <- c(curve_params, vt_params, restriction_levels, betas)

  restrictions_params <- cfg$restrictions_params
  restriction_st_params <- cfg$restriction_st_params
  param_scale <- abs(ini_params) / 10
  waves_list <- ranges_to_waves(cfg$waves_list)
  params_limits <- cfg$params_limits
  num_waves <- cfg$num_waves

  waves <- waves_1d_list(num_waves, waves_list)

  rho <- eval(parse(text = cfg$rho))
  infections <- cases / rho # Estimated time series of total infections

  # ==============================================================================

    est <- estimate.mle(
      ini_params = ini_params,
      params_limits = params_limits,
      restrictions = restrictions_params,
      restriction.starts = restriction_st_params,
      ranges = waves,
      rt_func = rt.func.num,
      silence.errors = silence.errors,
      fit.t.pred = fit.t.pred,
      param_scale = param_scale,
      num.iter = num.iter,
      cases = cases,
      scenario = NULL,
      H.E = NULL,
      H.W = NULL,
      adj.period = adj.period,
      population = population,
      rho = rho,
      serial_mean = serial_mean,
      serial_var = serial_var,
      lt = lt,
      window_size = window_size,
      hessian = T
    )

    print(est)

    a1 <- est$a1
    a2 <- est$a2
    a3 <- est$a3
    a4 <- est$a4
    nu <- est$nu
    vt <- c(1, est$vt_params_est)
    psi <- est$Psi
    betas <- est$betas


    # Result & Plot: Predictions are based on the function below

    r1 <- pred.curve(
      a1 = a1,
      a2 = a2,
      a3 = a3,
      a4 = a4,
      nu = nu,
      Psi = psi,
      betas = betas,
      use.actual.not.predicted = not.predict,
      restrictions = restrictions_params,
      restriction.starts = restriction_st_params,
      ranges = waves,
      variant.transm = vt,
      rt_func = rt.func.num,
      fit.t.pred = fit.t.pred,
      cases = cases,
      scenario = NULL,
      H.E = NULL,
      H.W = NULL,
      adj.period = adj.period,
      population = population,
      rho = rho,
      serial_mean = serial_mean,
      serial_var = serial_var,
      lt = lt,
      window_size = window_size
    )

    print(r1)
    Rt.pr <- r1$`Predicted Rt` # Predicted Rt
    yt.pr1 <- abs(r1$`Predicted Infections`)

    # Rt
    range.rt <- 1:length(lt)

    ### Plot for each Rt (1:4)
    range.points <- 1:(max(waves_list[[4]]) + 21) # Decide the range of points to plot

    # Rt
    plot_outputs(Time = Time,
                 cases = cases,
                 window_size = window_size,
                 serial_mean = serial_mean,
                 serial_var = serial_var,
                 predict.beyond = predict.beyond,
                 waves_list = waves_list,
                 num_waves = num_waves,
                 rt_func = rt.func.num,
                 curve = r1,
                 restrictions = restrictions_params,
                 restriction.starts = restriction_st_params
    )

    # Confidence bounds
    bounds <- ci.curve(fit = est,
                       restrictions = restrictions_params,
                       restriction.starts = restriction_st_params,
                       ranges = waves,
                       rt_func = rt.func.num,
                       fit.t.pred = fit.t.pred,
                       cases = cases,
                       scenario = NULL,
                       H.E = NULL,
                       H.W = NULL,
                       adj.period = adj.period,
                       population = population,
                       rho = rho,
                       serial_mean = serial_mean,
                       serial_var = serial_var,
                       lt = lt,
                       window_size = window_size)
    print(bounds)
}
