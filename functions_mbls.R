# Functions:

# custom function to get national level forecasts for a given target
# path: the folder containing all forecast files
# tg: the target
get_forecasts <- function(path, tg, season = NA){
  forecast_files0 <- list.files(path) # get fle names
  extr_dat <- function(str) strsplit(str, "201")[[1]][2]
  # re-order:
  order_forecast_files <- order(sapply(forecast_files0, extr_dat))
  forecast_files <- forecast_files0[order_forecast_files]

  issue_weeks <- sapply(forecast_files, extract_issue_week_from_filename)

  # load the forecasts from the first week to get dimensions:
  dat1 <- read.csv(paste0(path, "/", forecast_files[[1]]))
  colnames(dat1) <- tolower(colnames(dat1))
  dat_pw1 <- subset(dat1, target == tg &
                      location == "US National" &
                      type == "Bin")

  # now get all forecasts:
  all_forecasts <- matrix(ncol = nrow(dat_pw1), nrow = length(forecast_files))
  rownames(all_forecasts) <- paste0("iw", issue_weeks)
  colnames(all_forecasts) <- paste0("bin", as.numeric(as.character(dat_pw1$bin_start_incl)))

  for(i in seq_along(forecast_files)){
    dat <- read.csv(paste0(path, "/", forecast_files[[i]]))
    colnames(dat) <- tolower(colnames(dat))
    dat_pw <- subset(dat, target == tg &
                       location == "US National" &
                       type == "Bin")
    all_forecasts[i, ] <- dat_pw$value
  }
  return(list(forecast = all_forecasts,
              issue_weeks = issue_weeks,
              season = season,
              target = tg,
              optimized = FALSE))
}

# compute expected log score (uni-bin)
# pred: the prediction (numeric vector)
# true_distr: the true distribution (numeric vector)
e_ls <- function(pred, true_distr){
  score_values <- pred
  weighted_score_values <- true_distr*log(score_values)
  expected_score <- sum(weighted_score_values[true_distr > 0])
  return(expected_score)
}

# blur distribution, i.e. move from F to F_tilde
# distr: a distribution over ordered categories (numeric vector)
# d: degree of blurring (spread out over how many neigbouring bins?)
blur_one_forecast <- function(distr, d){
  lgt <- length(distr)
  distr_tilde <- numeric(lgt)
  # compute separately for each level:
  for(i in 1:lgt){
    subs_temp <- max(1, i - d):min(lgt, i + d)
    distr_tilde[i] <- sum(distr[subs_temp])/(2*d + 1)
  }
  return(distr_tilde)
}

# apply blurring to each row of a data frame
# forecast: data frame containing one distribution per row
# d: maximum acceptable difference
blur_forecasts <- function(forecast, d){
  ret <- list()
  ret$forecast <- t(apply(forecast$forecast, 1, blur_one_forecast, d = d))
  rownames(ret$forecast) <- rownames(forecast$forecast)
  colnames(ret$forecast) <- colnames(forecast$forecast)
  ret$issue_weeks <- forecast$issue_weeks
  ret$season = forecast$season
  ret$target <- forecast$target
  return(ret)
}

# compute expected multi-bin log score
# pred: the prediction (numeric vector)
# true_distr: the true distribution (numeric vector)
# d: maximum acceptable difference (number of bins on either side counted as correct)
e_mbls <- function(pred, true_distr, d = 1){
  score_values <- (2*d + 1)*blur_one_forecast(distr = pred, d = d) # compute mbls values
  weighted_score_values <- true_distr*log(score_values) # weigh with true probabilities
  expectation <- sum(weighted_score_values[true_distr > 0]) # compute expectation
  return(expectation)
}

# compute expected multi-bin log score for each row of a data frame
# pred_df: data frame containing one predictive distribution per row
# true_df: data frame containing the true distributions
# d: maximum acceptable difference
e_mbls_df <- function(pred_df, true_df, d){
  n_timepoints <- nrow(pred_df)
  res <- numeric(n_timepoints)
  for(i in 1:n_timepoints){
    res[i] <- e_mbls(pred = pred_df[i, ],
                     true_distr = true_df[i, ],
                     d = d)
  }
  return(res)
}

# anti-logit function (needed in numerical optimization)
antilogit <- function(p){
  sum_of_exp <- exp(1)/p[length(p)]
  exps <- p[-length(p)]*sum_of_exp
  log(exps)
}

# find forecast G which optimizes mbls under F
# original_forecast: the original forecast F (numeric vector)
# detailed: should the entire return object from optim be returned?
# d: maximum acceptable difference
# control_optim: optional control arguments passed to optim
optimize_one_forecast <- function(original_forecast, d = 1, detailed = TRUE, control_optim = NULL){

  original_forecast_long <- original_forecast
  relevant_range <- min(which(original_forecast > 1/(10*length(original_forecast)))):
    max(which(original_forecast > 1/(10*length(original_forecast))))
  original_forecast_short <- original_forecast[relevant_range]

  # function which will be minimized:
  # par: predictive distribution as vector on a multinomial anti-logit scale
  # note: length reduced by one as p_T follows from normalization
  fct_to_optimize <- function(par){
    # obtain predictive distribution on probability scale
    pred <- c(exp(par), 1)/sum(c(exp(par), 1))
    # return negative expected mbls
    -1*e_mbls(pred = pred, true_distr = original_forecast_short, d = d)
  }

  # use mixture of original forecast distribution and uniform as starting value:
  original_forecast_short_plus_eps <- 0.99*original_forecast_short +
    rep(0.01/length(original_forecast), length(original_forecast_short))
  initial <- antilogit(original_forecast_short_plus_eps)[-length(original_forecast_short_plus_eps)]
  initial <- pmax(-10, initial) # avoid too small values

  # run optimization
  opt <- optim(par = initial, fct_to_optimize, control = control_optim)
  # transform result back to probability scale:
  best_pred_short <- c(exp(opt$par), 1)/sum(c(exp(opt$par), 1))
  best_pred <- original_forecast
  best_pred[relevant_range] <- best_pred_short*sum(original_forecast_short) # normalization

  if(!detailed){
    return(best_pred)
  }else{
    return(list(forecast = best_pred,
                e_mbls = -opt$value,
                convergence = opt$convergence,
                opt = opt))
  }
}

# apply optimize_one_forecast to all rows of a data frame
# original_forecast: data frame containing one true distribution per row
# d: maximum acceptable difference
optimize_forecasts <- function(original_forecast, d){
  n_forecasts <- nrow(original_forecast$forecast)

  all_tuned_forecasts <- original_forecast$forecast*NA
  convergence <- e_mbls <- numeric(n_forecasts)
  pb <- txtProgressBar(min = 0, max = n_forecasts, style = 3)

  for(i in 1:n_forecasts){
    try({
      optimized_forecast_temp <- optimize_one_forecast(original_forecast$forecast[i, ],
                                                       d = d,
                                                       control_optim = list(maxit = 100000))
      all_tuned_forecasts[i, ] <- optimized_forecast_temp$forecast
      convergence[i] <- optimized_forecast_temp$convergence
      e_mbls[i] <- optimized_forecast_temp$e_mbls

      setTxtProgressBar(pb, i)
    })
  }
  return(list(forecast = all_tuned_forecasts,
              issue_weeks = original_forecast$issue_weeks,
              target = original_forecast$target,
              e_mbls = e_mbls,
              optimized = TRUE,
              convergence = convergence))
}

# helper function to extract the week number from a file name
extract_issue_week_from_filename <- function(filename){
  as.numeric(paste(strsplit(filename, split = "", fixed = TRUE)[[1]][3:4], collapse = ""))
}

# evaluate multi-bin log score for a forecast
# forecast: the forecast distribution, a numeric vector
# original_forecast
evaluate_one_forecast_mbls <- function(forecast, observed_bin, d){
  ind_observed_bin <- which(names(forecast) == paste0("bin", observed_bin))
  inds_acceptable_bins <- max(1, ind_observed_bin - d):
    min(length(forecast), ind_observed_bin + d)
  log(sum(forecast[inds_acceptable_bins]))
}

# apply evaluate_one_forecast_mbls to each row of a data.frame
evaluate_forecasts_mbls <- function(forecast, observed_bin, d){
  n_forecasts <- nrow(forecast$forecast)
  mbls <- numeric(n_forecasts)
  names(mbls) <- rownames(forecast$forecast)
  for(i in 1:n_forecasts){
    mbls[i] <- evaluate_one_forecast_mbls(forecast = forecast$forecast[i, ],
                                          observed_bin = observed_bin[i],
                                          d = d)
  }
  return(mbls)
}

# plotting function
plot_mbls <- function(timepoints, mbls_F, mbls_G, col = c("red", "blue"), legend = TRUE, pos.legend = "bottomright",...){
  plot(timepoints, mbls_F, type = "l", xlab = "week in which forecast is issued", ylab = "MBlogS", axes = FALSE,
       col = col[1], ylim = c(-2.5, 0), xlim = c(0, 33),...)
  axis(1, at = c(1, 7, 14, 21, 28))
  # axis(3, at = c(1, 7, 14, 21, 28), labels = rownames(mbls_F)[c(1, 7, 14, 21, 28)])
  # mtext(3, line = 3, text = "Calendar week")
  axis(2)
  box()
  lines(timepoints, mbls_G, col = col[2])
  if(legend){
    legend(pos.legend, legend = c(paste0("F; mean MBlogS: ", round(mean(mbls_F), 3)),
                                  paste("G; mean MBlogS: ", round(mean(mbls_G), 3))),
           col = c(col), lty = 1, bty = "n")
  }
}
