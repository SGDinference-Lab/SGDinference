#' Averaged SGD and its Inference via Random Scaling
#'
#' Compute the averaged SGD estimator and the confidence intervals via random scaling method.
#'
#' @param formula formula. The response is on the left of a ~ operator. The terms are on the right of a ~ operator, separated by a + operator.
#' @param data an optional data frame containing variables in the model. 
#' @param gamma_0 numeric
#' @param alpha numeric
#' @param burn numeric
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
#' my.dat = data.frame(y=y, x=x)
#' sgd.out = sgd_lm(y~., data=my.dat)


# Todo list
# (1) "rss" subset inference for linear regression

sgd_lm = function(formula, data, gamma_0=1, alpha=0.667, burn=1, 
                bt_start = NULL,  
                studentize = TRUE, intercept = TRUE,
                inference = "rs"
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
  # Linear (Mean) Regression 
  #----------------------------------------------
  out = sgd_lm_cpp(x, y, burn, gamma_0, alpha, bt_start)
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
result.out$coefficient = beta_hat
result.out$call = cl
result.out$terms <- mt
result.out$var <- NULL

result.out$ci.lower = NULL
result.out$ci.upper = NULL
  
return(result.out)

}


# beta_hat_path and V_hat_path are not completed in this code.
