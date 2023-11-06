#' Averaged S-subGD Estimator in Linear Quantile Regression
#'
#' Compute the averaged S-subGD (stochastic subgradient) estimator for the coefficients in linear quantile regression.
#'
#' @param formula formula. The response is on the left of a ~ operator. The terms are on the right of a ~ operator, separated by a + operator.
#' @param data an optional data frame containing variables in the model. 
#' @param gamma_0 numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is NULL and it is determined by the adaptive method in Lee et al. (2023).
#' @param alpha numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 0.501.
#' @param burn numeric. A tuning parameter for "burn-in" observations. 
#'    We burn-in up to (burn-1) observations and use observations from (burn) for estimation. Default is 1, i.e. no burn-in. 
#' @param bt_start numeric. (p x 1) vector, excluding the intercept term. User-provided starting value. Default is NULL. Then, it is estimated by conquer.
#' @param qt numeric. Quantile. Default is 0.5. 
#' @param studentize logical. Studentize regressors. Default is TRUE.
#' @param no_studentize numeric. The number of observations to compute the mean and std error for studentization. Default is 100. 
#' @param intercept logical. Use the intercept term for regressors. Default is TRUE. 
#'    If this option is TRUE, the first element of the parameter vector is the intercept term.
#' @param path logical. The whole path of estimation results is out. Default is FALSE.
#' @param path_index numeric. A vector of indices to print out the path. Default is 1.
#'
#' @return
#' An object of class \code{"sgdi"}, which is a list containing the following
#' \describe{
#' \item{\code{coefficients}}{a vector of estimated parameter values}
#' \item{\code{path_coefficients}}{The path of coefficients.}
#' }
#' @note{The dimension of \code{coefficients} is (p+1) if \code{intercept}=TRUE or p otherwise.}
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
#' sgd.out = sgd_qr(y~., data=my.dat)

sgd_qr = function(formula, 
                  data, 
                  gamma_0=NULL, 
                  alpha=0.501, 
                  burn=1, 
                  bt_start = NULL, 
                  qt=0.5,
                  studentize = TRUE, 
                  no_studentize = 100L,
                  intercept = TRUE,
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

  # Get the dimension of x and the sample size: p and n
  p = ncol(as.matrix(x))
  n = length(y)

  # Select gamma_0 by the data adaptive method
  if (is.null(gamma_0)){
    sig_hat = sd(y)
    gamma_0 = (dnorm(qnorm(qt))/sqrt(qt*(1-qt))) / sig_hat
  }
  
  # Initialize the bt_t
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

  #----------------------------------------------
  # Quantile Regression
  #----------------------------------------------

    out = sgd_qr_cpp(x, y, burn, gamma_0, alpha, bt_start=bt_start, tau=qt, x_mean=x_mean_in, x_sd=x_sd_in, path=path, path_index=path_index)
    beta_hat = out$beta_hat

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
  }


    
#--------------------------------------------
# out: list of all outputs
#--------------------------------------------
result.out = list()
class(result.out) = "sgdi"
result.out$coefficients = beta_hat
result.out$intercept = intercept
result.out$level = 0
result.out$call = cl
result.out$terms <- mt
result.out$V <- NULL
result.out$ci.lower = NULL
result.out$ci.upper = NULL
result.out$gamma_0 = gamma_0

if (path){
  if (studentize){
    result.out$path_coefficients = (out$beta_hat_path) %*% t(rescale_matrix[path_index, path_index])
  } else {
    result.out$path_coefficients = out$beta_hat_path
  }
}
return(result.out)

}
