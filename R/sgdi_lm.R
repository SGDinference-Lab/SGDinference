#' Averaged SGD and its Inference via Random Scaling
#'
#' Compute the averaged SGD estimator and conduct inference via random scaling method.
#'
#' @param formula formula. The response is on the left of a ~ operator. The terms are on the right of a ~ operator, separated by a + operator.
#' @param data an optional data frame containing variables in the model. 
#' @param gamma_0 numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 1.
#' @param alpha numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 0.667.
#' @param burn numeric. A tuning parameter for "burn-in" observations. We burn-in up to (burn-1) observations and use observations from (burn) for estimation. Default is 1, i.e. no burn-in. 
#' @param inference character. Specifying the inference method. Default is "rs" (random scaling matrix for joint inference using all the parameters). 
#'    "rss" is for ransom scaling subset inference. This option requires that "rss_indx" should be provided.
#'    "rsd" is for the diagonal elements of the random scaling matrix, excluding one for the intercept term.  
#' @param bt_start numeric. (p x 1) vector. User-provided starting value Default is NULL.
#' @param studentize logical. Studentize regressors. Default is TRUE
#' @param no_studentize numeric. The number of observations to compute the mean and std error for studentization. Default is 100. 
#' @param intercept logical. Use the intercept term for regressors. Default is TRUE. 
#'    If this option is TRUE, the first element of the parameter vector is the intercept term.
#' @param rss_idx numeric. Index of x for random scaling subset inference. Default is 1, the first regressor of x. 
#'    For example, if we want to focus on the 1st and 3rd covariates of x, then set it to be c(1,3).
#' @param level numeric. The confidence level required. Default is 0.95. Can choose 0.90 and 0.80.
#'
#' @return
#' An object of class \code{"sgdi"}, which is a list containing the following
#' \describe{
#' \item{\code{coefficient}}{A (p + 1)-vector of estimated parameter values including the intercept.}
#' \item{\code{var}}{A (p+1)x (p+1) variance-covariance matrix of \code{coefficient}}
#' \item{\code{ci.lower}}{The lower part of the 95\% confidence interval}
#' \item{\code{ci.upper}}{The upper part of the 95\% confidence interval}
#' }
#' @export
#'
#' @examples
#' n = 1e05
#' p = 5
#' bt0 = rep(5,p)
#' x = matrix(rnorm(n*(p-1)), n, (p-1))
#' y = cbind(1,x) %*% bt0 + rnorm(n)
#' my.dat = data.frame(y=y, x=x)
#' sgdi.out = sgdi_lm(y~., data=my.dat)

sgdi_lm = function(formula, 
                  data, 
                  gamma_0=1, 
                  alpha=0.667, 
                  burn=1, 
                  inference="rs",
                  bt_start = NULL,  
                  studentize = TRUE, 
                  no_studentize = 100L,
                  intercept = TRUE,
                  rss_idx = c(1), 
                  level = 0.95){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)  # model.frame returns
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- as.matrix(model.matrix(mt, mf, drop=F)[,-1])
  if (inference == "rss"){
    if (0 %in% rss_idx ){
      stop("rss_idx includes 0 (the intercept term), where it should be bigger than 1.")
    }
    rss_idx_r = rss_idx + 1   # Index starting with 1 in R and with 0 in C    
  }

  
  if (studentize){
    # Compute column means and standard errors and save them for later reconversion
    if (no_studentize > length(y)) {
      cat("Warning: no_studentize is bigger than the sample size. no_studentize is set to be the sample size. \n")
      no_studentize = length(y)
    }
    x_mean = apply(x[1:no_studentize, , drop=F], 2, mean)
    x_sd = apply(x[1:no_studentize, , drop=F], 2, sd)
    if (intercept){
      x = cbind(1, x)
      x_mean_in = c(0,x_mean)
      x_sd_in = c(1,x_sd)
    } else {
      x_mean_in = x_mean
      x_sd_in = x_sd
    }
  } else {
    if (intercept){
      x = cbind(1, x)
    }  
    # They are irrelevant in computation but should be initialized. We set some extreme numbers.
    x_mean_in = -1e06
    x_sd_in = -1.0
  }

  # Attach a vector of 1's for an intercept term

  # Get the dimension of x and the sample size: p and n
  p = ncol(as.matrix(x))
  n = length(y)

  # Initialize the bt_t, A_t, b_t, c_t
  if (is.null(bt_start)){
    bt_t = bar_bt_t = bt_start = matrix(0, nrow=p, ncol=1)
  } else {
    bt_t = bar_bt_t = matrix(bt_start, nrow=p, ncol=1)
  }
  A_t = matrix(0, p, p)
  b_t = matrix(0, p, 1)
  c_t = 0
  V_t = NULL


  #----------------------------------------------
  # Linear (Mean) Regression 
  #----------------------------------------------
  out = sgdi_lm_cpp(x, y, burn, gamma_0, alpha, bt_start=bt_t, inference=inference, rss_idx=rss_idx, x_mean=x_mean_in, x_sd=x_sd_in)
  beta_hat = out$beta_hat
  if (inference == "rs"){
    V_out = out$V_hat
  } else if (inference == "rss"){
    V_out = out$V_hat_sub
  } else if (inference == "rsd"){
    V_out = out$V_hat_diag
  }

  # Re-scale parameters to reflect the studentization
  if (studentize){
    if (length(x_sd)>1){
      rescale_matrix = diag(1/x_sd)
    } else {
      rescale_matrix = 1/x_sd
    }
    if (intercept){
      # Redefine the rescale_matrix including the intercept term
      rescale_matrix = rbind(c(1,-(x_mean/x_sd)), cbind(0, rescale_matrix))
    }
    # Re-scale the parameters
    beta_hat = rescale_matrix %*% beta_hat
    # Re-scale the variance
    if (inference == "rs"){
      V_out = rescale_matrix %*% V_out %*% t(rescale_matrix)  
    } else if (inference == "rss"){
      V_out = rescale_matrix[rss_idx_r, rss_idx_r] %*% V_out %*% t(rescale_matrix[rss_idx_r, rss_idx_r])  
    } else if (inference == "rsd"){
      V_out = diag(rescale_matrix) * V_out * diag(rescale_matrix)
      # With studentization, we cannot compute the variance of an intercept, which requires the whole V_hat.
      V_out[1] = NA 
    }
  }

  #--------------------------------------------
  # out: list of all outputs
  #--------------------------------------------
  result.out = list()
  class(result.out) = "sgdi"
  result.out$coefficients = beta_hat
  result.out$call = cl
  result.out$terms <- mt
  result.out$V <- V_out
    
  if (level == 0.95) {
    critical.value = 6.747       # From Abadir and Paruolo (1997) Table 1. 97.5%  
  } else if (level == 0.9) {
    critical.value = 5.323       # From Abadir and Paruolo (1997) Table 1. 95.0%  
  } else if (level == 0.8) {
    critical.value = 3.875       # From Abadir and Paruolo (1997) Table 1. 90.0%  
  } else {
    critical.value = 6.747
    cat("Confidence level should be chosen from 0.95, 0.90, and 0.80. \n")
    cat("We report the default level 0.95. \n")
  }
  
  if (inference == "rs"){
    ci.lower = beta_hat - critical.value * sqrt(diag(V_out)/n)
    ci.upper = beta_hat + critical.value * sqrt(diag(V_out)/n) 
  } else if (inference == "rss"){
    ci.lower = ci.upper = matrix(NA, p, 1)
    ci.lower[rss_idx_r] = beta_hat[rss_idx_r] - critical.value * sqrt(diag(V_out)/n)
    ci.upper[rss_idx_r] = beta_hat[rss_idx_r] + critical.value * sqrt(diag(V_out)/n) 
  } else if (inference == "rsd"){
    ci.lower = beta_hat - critical.value * sqrt(V_out/n)
    ci.upper = beta_hat + critical.value * sqrt(V_out/n) 
  }
  
  result.out$ci.lower = ci.lower
  result.out$ci.upper = ci.upper
  
  result.out$intercept = intercept
  result.out$inference = inference
  result.out$level = level
  
  if (inference == "rss"){
    result.out$rss_idx_r = rss_idx_r
  }
  
  return(result.out)
  
}
  
  
  