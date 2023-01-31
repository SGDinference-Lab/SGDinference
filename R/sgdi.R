#' Averaged SGD and its Inference via Random Scaling
#'
#' Compute the averaged SGD estimator and the variance-covariance matrix via random scaling method.
#'
#' @param formula formula. The response is on the left of a ~ operator. The terms are on the right of a ~ operator, separated by a + operator.
#' @param data an optional data frame containing variables in the model. 
#' @param gamma_0 numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 1.
#' @param alpha numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 0.667.
#' @param burn numeric. A tuning parameter for "burn-in" observations. We burn-in up to (burn-1) observations and use observations from (burn) for estimation. Default is 1, i.e. no burn-in. 
#' @param model character specifying the model to be used: \code{"lm"} (linear
#'   mean regression model), \code{"qr"} (quantile regression)
#' @param inference character. Specifying the inference method. Default is "rs" (random scaling). "rss" is for ransom scaling subset inference. Then, "rss_indx" should be provided. 
#' @param bt_start numeric. (p x 1) vector. User-provided starting value Default is NULL.
#' @param path_output numeric specifying the sequence that print out the output paths
#' @param qt numeric. Quantile. Default is 0.5. 
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
#' my.dat = data.frame(y=y, x=x)
#' out = sgdi(y~., data=my.dat)


sgdi = function(formula, data, gamma_0=1, alpha=0.667, burn=1, model="lm", inference="rs",
                     bt_start = NULL, path_output = NULL, qt=0.5,
                     studentize = TRUE, intercept = TRUE, rss_idx = c(1)
){
  V_hat_sub = 0
  cl <- match.call()
  
  #----------------------------------------------
  # Linear (Mean) Regression 
  #----------------------------------------------
  if (model=="lm"){
    out = sgdi_lm(formula, data, gamma_0, alpha, burn, inference,
                   bt_start, path_output, 
                   studentize, intercept
    )
  }
  
  
  #----------------------------------------------
  # Quantile Regression
  #----------------------------------------------
  if (model=="qr"){
    out = sgdi_qr(formula, data, gamma_0, alpha, burn, inference,
                  bt_start, path_output, qt,
                  studentize, intercept
    )
  }
  
  out$call = cl
  
  if ( is.null(path_output)) {
    return(out)
  } else {
    NULL
    # The path return was note written yet.
    # return(list(beta_hat = beta_hat, V_hat = V_hat, beta_hat_path = beta_hat_path, V_hat_path = V_hat_path))
  }
  
}



