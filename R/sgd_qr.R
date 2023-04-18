#' Averaged S-subGD Estimator in Linear Quantile Regression
#'
#' Compute the averaged S-subGD (stochastic subgradient) estimator for the coefficients in linear quantile regression.
#'
#' @param formula formula. The response is on the left of a ~ operator. The terms are on the right of a ~ operator, separated by a + operator.
#' @param data an optional data frame containing variables in the model. 
#' @param gamma_0 numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 1.
#' @param alpha numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 0.667.
#' @param burn numeric. A tuning parameter for "burn-in" observations. 
#'    We burn-in up to (burn-1) observations and use observations from (burn) for estimation. Default is 1, i.e. no burn-in. 
#' @param bt_start numeric. (p x 1) vector, excluding the intercept term. User-provided starting value. Default is NULL.
#' @param qt numeric. Quantile. Default is 0.5. 
#' @param studentize logical. Studentize regressors. Default is TRUE.
#' @param intercept logical. Use the intercept term for regressors. Default is TRUE. 
#'    If this option is TRUE, the first element of the parameter vector is the intercept term.
#'
#' @return
#' An object of class \code{"sgdi"}, which is a list containing the following
#' \describe{
#' \item{\code{coefficients}}{a vector of estimated parameter values}
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

sgd_qr = function(formula, data, gamma_0=1, alpha=0.667, burn=1, 
                bt_start = NULL, qt=0.5,
                studentize = TRUE, intercept = TRUE
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

  #----------------------------------------------
  # Quantile Regression
  #----------------------------------------------

    out = sgd_qr_cpp(x, y, burn, gamma_0, alpha, bt_start=bt_t, tau=qt)
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
result.out$call = cl
result.out$terms <- mt
result.out$V <- NULL
result.out$ci.lower = NULL
result.out$ci.upper = NULL
  
return(result.out)

}
