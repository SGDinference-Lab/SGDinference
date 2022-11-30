#' Averaged SGD and its Inference via Random Scaling
#'
#' Compute the averaged SGD estimator and the variance-covariance matrix of linear mean regression via random scaling method.
#'
#' @param x numeric. (n x p) matrix of regressors. Should not include 1 (the intercept)
#' @param y numeric. (n x 1) vector of a dependent variable. 
#' @param gamma_0 numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 1.
#' @param alpha numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 0.667.
#' @param burn numeric. A tuning parameter for "burn-in" observations. We burn-in up to (burn-1) observations and use observations from (burn) for estimation. Default is 1, i.e. no burn-in. 
#' @param inference character. Specifying the inference method. Default is "rs" (random scaling). "rss" is for ransom scaling subset inference. Then, "rss_indx" should be provided. 
#' @param bt_start numeric. (p x 1) vector. User-provided starting value Default is NULL.
#' @param path_output numeric specifying the sequence that print out the output paths
#' @param studentize logical. Studentize regressors. Default is TRUE
#' @param intercept logical. Use the intercept term for regressors. Default is TRUE
#' @param rss_idx numeric. Index of x for random scaling subset inference. Default is 1, the first regressor of x. For example, if we want to infer the 1st, 3rd covariate of x, then set it to be c(1,3).
#'
#' @return
#' An object of class \code{"sgdi"}, which is a list containing the following
#' \describe{
#' \item{\code{beta_hat}}{A (p + 1)-vector of estimated parameter values including the intercept.}
#' \item{\code{V_hat}}{A (p+1)x (p+1) variance-covariance matrix of \code{beta_hat}}
#' \item{\code{V_hat_sub}}{A variance-covariance sub-matrix of \code{beta_hat}. If the subset size is not provided, it returns 0.}
#' }
#' @export
#'
#' @examples
#' n = 1e05
#' p = 5
#' bt0 = rep(5,p)
#' x = matrix(rnorm(n*(p-1)), n, (p-1))
#' y = cbind(1,x) %*% bt0 + rnorm(n)
#' sgdi.out = sgdi_lm(x,y)
#' out = sgdi_lm(x,y)

# Todo list
# (1) "rss" subset inference for linear regression
# (2) path_output

sgdi_lm = function(x, y, gamma_0=1, alpha=0.667, burn=1, inference="rs",
                bt_start = NULL, path_output = NULL, 
                studentize = TRUE, intercept = TRUE
                ){
  x = as.matrix(x)
  
  if (studentize){
    # Compute column means and standard errors and save them for later reconversion
    x_mean = apply(x, 2, mean)
    x_sd = apply(x, 2, sd)

    # Studentize each column if x
    x = apply(x, 2, function(.) (.-mean(.))/sd(.) )
  }

  # Attach a vector of 1's for an intercept term
  if (intercept){
    x = cbind(1, x)
  }

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

  n_path = length(path_output)
  cnt_path = 1
  beta_hat_path = matrix(NA, p, n_path)
  V_hat_path = array(NA, dim = c(p, p, n_path))

  #----------------------------------------------
  # Linear (Mean) Regression 
  #----------------------------------------------
  out = sgdi_lm_cpp(x, y, burn, gamma_0, alpha, bt_start, "rs")
  beta_hat = out$beta_hat
  V_hat = out$V_hat

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
    V_hat = rescale_matrix %*% V_hat %*% t(rescale_matrix)
  }



  if ( is.null(path_output)) {
    return(list(beta_hat=beta_hat, V_hat = V_hat))
  } else {
    return(list(beta_hat = beta_hat, V_hat = V_hat, beta_hat_path = beta_hat_path, V_hat_path = V_hat_path))
  }

}


# beta_hat_path and V_hat_path are not completed in this code.
