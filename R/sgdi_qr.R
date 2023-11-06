#' Averaged S-subGD and its Inference via Random Scaling in Linear Quantile Regression
#'
#' Compute the averaged S-subGD (stochastic subgradient) estimator for the coefficients in linear quantile regression and conduct inference via random scaling method.
#'
#' @param formula formula. The response is on the left of a ~ operator. The terms are on the right of a ~ operator, separated by a + operator.
#' @param data an optional data frame containing variables in the model. 
#' @param gamma_0 numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is NULL and it is determined by the adaptive method in Chet et al. (2023).
#' @param alpha numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 0.501.
#' @param burn numeric. A tuning parameter for "burn-in" observations. 
#'    We burn-in up to (burn-1) observations and use observations from (burn) for estimation. Default is 1, i.e. no burn-in. 
#' @param inference character. Specifying the inference method. Default is "rs" (random scaling matrix for joint inference using all the parameters). 
#'    "rss" is for ransom scaling subset inference. This option requires that "rss_indx" should be provided.
#'    "rsd" is for the diagonal elements of the random scaling matrix, excluding one for the intercept term.  
#' @param bt_start numeric. (p x 1) vector, excluding the intercept term. User-provided starting value. Default is NULL. Then, it is estimated by conquer.
#' @param qt numeric. Quantile. Default is 0.5. 
#' @param studentize logical. Studentize regressors. Default is TRUE.
#' @param no_studentize numeric. The number of observations to compute the mean and std error for studentization. Default is 100. 
#' @param intercept logical. Use the intercept term for regressors. Default is TRUE. 
#'    If this option is TRUE, the first element of the parameter vector is the intercept term.
#' @param rss_idx numeric. Index of x for random scaling subset inference. Default is 1, the first regressor of x. 
#'    For example, if we want to focus on the 1st and 3rd covariates of x, then set it to be c(1,3).
#' @param level numeric. The confidence level required. Default is 0.95. Can choose 0.90 and 0.80. 
#' @param path logical. The whole path of estimation results is out. Default is FALSE.
#' @param path_index numeric. A vector of indices to print out the path. Default is 1.
#'
#' @return
#' An object of class \code{"sgdi"}, which is a list containing the following
#' \describe{
#' \item{\code{coefficients}}{a vector of estimated parameter values}
#' \item{\code{V}}{a random scaling matrix depending on the inference method}
#' \item{\code{ci.lower}}{a vector of lower confidence limits}
#' \item{\code{ci.upper}}{a vector of upper confidence limits}
#' \item{\code{inference}}{character that specifies the inference method}
#' }
#' @note{The dimension of \code{coefficients} is (p+1) if \code{intercept}=TRUE or p otherwise.
#' The random scaling matrix \code{V} is a full matrix if "rs" is chosen;
#' it is a scalar or smaller matrix, depending on the specification of "rss_indx" if "rss" is selected;
#' it is a vector of diagonal elements of the full matrix if "rsd" is selected. 
#' In this case, the first element is missing if the intercept is included.
#' The confidence intervals may contain NA under "rss" and "rsd".}
#' 
#' @export
#'
#' @examples
#' n = 1e05
#' p = 5
#' bt0 = rep(5,p)
#' x = matrix(rnorm(n*(p-1)), n, (p-1))
#' y = cbind(1,x) %*% bt0 + rnorm(n)
#' my.dat = data.frame(y=y, x=x)
#' sgdi.out = sgdi_qr(y~., data=my.dat)

sgdi_qr = function(formula, 
                  data, 
                  gamma_0=NULL, 
                  alpha=0.501, # Might need to fix 0.501
                  burn=1, 
                  inference="rs",
                  bt_start = NULL, 
                  qt=0.5,
                  studentize = TRUE, 
                  no_studentize = 100L,
                  intercept = TRUE,
                  rss_idx = c(1), 
                  level = 0.95,
                  path = FALSE,
                  path_index = c(1)){
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
  }
  # Attach a vector of 1's for an intercept term
  
  # Get the dimension of x and the sample size: p and n
  p = ncol(as.matrix(x))
  n = length(y)
  
  # Select gamma_0 by the data adaptive method
  if (is.null(gamma_0)){
    sig_hat = sd(y)
    gamma_0 = (dnorm(qnorm(qt))/sqrt(qt*(1-qt))) / sig_hat
  }
  
  # Initialize the bt_t, A_t, b_t, c_t
  if (is.null(bt_start)){
    # bt_t = bar_bt_t = bt_start = matrix(0, nrow=p, ncol=1)
    n_s = floor(max(c(n*0.01,p*10)))
    subsample_index = sample(n, n_s)
    if (intercept) {
      bt_start = conquer::conquer(x[subsample_index,-1, drop=F], y[subsample_index], tau = qt)$coeff 
    } else {
      bt_start = conquer::conquer(x[subsample_index, , drop=F], y[subsample_index], tau = qt)$coeff[-1]
    }
  } 
  A_t = matrix(0, p, p)
  b_t = matrix(0, p, 1)
  c_t = 0
  V_t = NULL
  
  #----------------------------------------------
  # Quantile Regression
  #----------------------------------------------
  out = sgdi_qr_cpp(x, y, burn, gamma_0, alpha, bt_start=bt_start, inference=inference, tau=qt, rss_idx=rss_idx, x_mean=x_mean_in, x_sd=x_sd_in, path=path, path_index=path_index)
  beta_hat = out$beta_hat
  if (inference == "rs"){
    V_out = out$V_hat
  } else if (inference == "rss"){
    V_out = out$V_hat_sub
  } else if (inference == "rsd"){
    V_out = out$V_hat_diag
  }
  
  if (studentize){
    # Re-scale the parameters
    beta_hat = rescale_matrix %*% beta_hat
    # Re-scale the variance
    if (inference == "rs"){
      V_out = rescale_matrix %*% V_out %*% t(rescale_matrix)  
    } else if (inference == "rss"){
      V_out = rescale_matrix[rss_idx_r, rss_idx_r] %*% V_out %*% t(rescale_matrix[rss_idx_r, rss_idx_r])  
    } else if (inference == "rsd"){
      V_out = diag(rescale_matrix) * V_out * diag(rescale_matrix)
      # If both "intercept" and "studentize" options are selected, 
      # we cannot compute the first diagonal element of the random scaling matrix because it requires the whole V_hat.
      if (intercept){
        V_out[1] = NA 
      }
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
  
  if (path){
    result.out$beta_hat_path = out$beta_hat_path
  }
  
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
  result.out$gamma_0 = gamma_0
  
  if (inference == "rss"){
    result.out$rss_idx_r = rss_idx_r
  }
  
  return(result.out)
  
}
