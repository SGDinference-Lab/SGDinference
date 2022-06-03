#' Averaged SGD and its Inference via Random Scaling
#'
#' Compute the averaged SGD estimator and the confidence intervals via random scaling method.
#'
#' @param x numeric. (n x p) matrix of regressors. Should not include 1 (the intercept)
#' @param y numeric
#' @param gamma_0 numeric
#' @param alpha numeric
#' @param burn numeric
#' @param path_output numeric specifying the sequence that print out the output paths
#' @param bt_start numeric
#' @param studentize logical. Studentize regressors. Default is TRUE
#' @param intercept logical. Use the intercept term for regressors. Default is TRUE
#' @param qt numeric. Quantile. Default is 0.5. 
#' @param h numeric. Bandwidth for smoothed QR estimator. Default is 0 (automatic choice)
#'
#' @return
#' #' An object of class \code{"sgdi"}, which is a list containing the following
#' components:
#'
#' @export
#'



sgd_sqr = function(x, y, gamma_0=1, alpha=0.667, burn=1, 
                bt_start = NULL, path_output = NULL, qt=0.5,
                studentize = TRUE, intercept = TRUE, h=0
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

  n_path = length(path_output)
  cnt_path = 1
  beta_hat_path = matrix(NA, p, n_path)

  #----------------------------------------------
  # Quantile Regression
  #----------------------------------------------

    out = sgd_sqr_cpp(x, y, burn, gamma_0, alpha, bt_start=bt_t, tau=qt, h=h)
    beta_hat = out$beta_hat

  # Re-scale parameters to reflect the studentization
  if (studentize){
    rescale_matrix = diag(1/x_sd)
    if (intercept){
      # Redefine the rescale_matrix including the intercept term
      rescale_matrix = rbind(c(1,-(x_mean/x_sd)), cbind(0, rescale_matrix))
    }
    # Re-scale the parameters
    beta_hat = rescale_matrix %*% beta_hat
  }



  if ( is.null(path_output)) {
    return(list(beta_hat=beta_hat))
  } else {
    return(list(beta_hat = beta_hat, beta_hat_path = beta_hat_path))
  }

}


# beta_hat_path and V_hat_path are not completed in this code.
