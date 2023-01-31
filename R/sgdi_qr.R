#' Inference of quantile regression with SGD and random scaling
#'
#' Compute the averaged SGD estimator of quantile regression and the confidence intervals of quantile regressionvia random scaling method.
#'
#' @param formula formula. The response is on the left of a ~ operator. The terms are on the right of a ~ operator, separated by a + operator.
#' @param data an optional data frame containing variables in the model. 
#' @param gamma_0 numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 1.
#' @param alpha numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 0.667.
#' @param burn numeric. A tuning parameter for "burn-in" observations. We burn-in up to (burn-1) observations and use observations from (burn) for estimation. Default is 1, i.e. no burn-in. 
#' @param inference character. Specifying the inference method. Default is "rs" (random scaling). "rss" is for ransom scaling subset inference. Then, "rss_indx" should be provided. 
#' @param bt_start numeric. (p x 1) vector. User-provided starting value Default is NULL.
#' @param path_output numeric specifying the sequence that print out the output paths
#' @param qt numeric. Quantile. Default is 0.5. 
#' @param studentize logical. Studentize regressors. Default is TRUE
#' @param intercept logical. Use the intercept term for regressors. Default is TRUE
#' @param rss_idx numeric. Index of x for random scaling subset inference. Default is 1, the first regressor of x. For example, if we want to infer the 1st, 3rd covariate of x, then set it to be c(1,3).
#'
#' @return
#' An object of class \code{"sgdi_qr"}, which is a list containing the following
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
#' sgdi.out = sgdi_qr(y~., data=my.dat)

# Todo list
# (1) path_output

sgdi_qr = function(formula, data, gamma_0=1, alpha=0.667, burn=1, inference="rs",
                bt_start = NULL, path_output = NULL, qt=0.5,
                studentize = TRUE, intercept = TRUE,
                rss_idx = c(1)
                ){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)  # model.frame returns
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)[,-1]
  
  
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

    out = sgdi_qr_cpp(x, y, burn, gamma_0, alpha, bt_start=bt_t, inference=inference, tau=qt, rss_idx=rss_idx)
    beta_hat = out$beta_hat
    V_hat = out$V_hat
    V_hat_sub = out$V_hat_sub



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
    V_hat_sub = rescale_matrix[rss_idx+1,rss_idx+1] %*% V_hat_sub %*% t(rescale_matrix[rss_idx+1,rss_idx+1])
  }



    #--------------------------------------------
    # out: list of all outputs
    #--------------------------------------------
    result.out = list()
    class(result.out) = "sgdi"
    result.out$coefficient = beta_hat
    result.out$call = cl
    result.out$terms <- mt
    result.out$var <- V_hat
    
    critical.value = 6.747       # From Abadir and Paruolo (1997) Table 1. 97.5%
    ci.lower = beta_hat - critical.value * sqrt(diag(V_hat)/n)
    ci.upper = beta_hat + critical.value * sqrt(diag(V_hat)/n) 
    result.out$ci.lower = ci.lower
    result.out$ci.upper = ci.upper
    
    if ( is.null(path_output)) {
      return(result.out)
    } else {
      return(list(beta_hat = beta_hat, V_hat = V_hat, beta_hat_path = beta_hat_path, V_hat_path = V_hat_path))
    }

}


# beta_hat_path and V_hat_path are not completed in this code.
