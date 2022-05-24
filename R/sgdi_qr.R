#' Averaged SGD and its Inference via Random Scaling
#'
#' Compute the averaged SGD estimator and the confidence intervals via random scaling method.
#'
#' @param x numeric. (n x p) matrix of regressors. Should not include 1 (the intercept)
#' @param y numeric
#' @param gamma_0 numeric
#' @param alpha numeric
#' @param burn numeric
#' @param model character specifying the model to be used: \code{"lm"} (linear
#'   model)
#' @param z numeric. (n x q) matrix of instruments used for \code{"tsls"} (Two-Stage Least Squares)
#' @param inference character specifying the inference method. Default is "rs" (random scaling)
#' @param path_output numeric specifying the sequence that print out the output paths
#' @param bar_Pi_s temporary. Will delete this later. (Population Pi_star)
#' @param bt_start numeric
#' @param studentize logical. Studentize regressors. Default is TRUE
#' @param intercept logical. Use the intercept term for regressors. Default is TRUE
#'
#' @return
#' #' An object of class \code{"sgdi"}, which is a list containing the following
#' components:
#'
#' @export
#'
#' @examples
#' n = 1e05
#' p = 5
#' bt0 = rep(5,p)
#' x = matrix(rnorm(n*(p-1)), n, (p-1))
#' y = cbind(1,x) %*% bt0 + rnorm(n)
#' sgdi.out = sgdi(x,y)


sgdi_qr = function(x, y, gamma_0=1, alpha=0.667, burn=1, inference="rs",
                bt_start = NULL, path_output = NULL, qt=0.5,
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
  # Quantile Regression
  #----------------------------------------------

    out = sgdi_qr_cpp(x, y, burn, gamma_0, alpha, bt_start=bt_t, "rs", tau=qt)
    beta_hat = out$beta_hat
    V_hat = out$V_hat



  # Re-scale parameters to reflect the studentization
  if (studentize){
    rescale_matrix = diag(1/x_sd)
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
