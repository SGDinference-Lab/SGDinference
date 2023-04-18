#' Averaged SGD (or S-subGD) and its Inference via Random Scaling in Linear Mean or Quantile Regression Models
#'
#' Compute the averaged SGD (or S-subGD) estimator and conduct inference via random scaling method.
#'
#'
#' @param formula formula. The response is on the left of a ~ operator. The terms are on the right of a ~ operator, separated by a + operator.
#' @param data an optional data frame containing variables in the model. 
#' @param gamma_0 numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 1.
#' @param alpha numeric. A tuning parameter for the learning rate (gamma_0 x t ^ alpha). Default is 0.667.
#' @param burn numeric. A tuning parameter for "burn-in" observations. 
#'    We burn-in up to (burn-1) observations and use observations from (burn) for estimation. Default is 1, i.e. no burn-in. 
#' @param model character specifying the model to be used: 
#'  \code{"lm"} (linear mean regression model), 
#'  \code{"qr"} (linear quantile regression model)
#' @param inference character. Specifying the inference method. Default is "rs" (random scaling matrix for joint inference using all the parameters). 
#'    "rss" is for ransom scaling subset inference. This option requires that "rss_indx" should be provided.
#'    "rsd" is for the diagonal elements of the random scaling matrix, excluding one for the intercept term.  
#' @param bt_start numeric. (p x 1) vector, excluding the intercept term. User-provided starting value. Default is NULL.
#' @param qt numeric. Quantile. Default is 0.5. This option is only relevant for quantile regression.
#' @param studentize logical. Studentize regressors. Default is TRUE.
#' @param intercept logical. Use the intercept term for regressors. Default is TRUE. 
#'    If this option is TRUE, the first element of the parameter vector is the intercept term.
#' @param rss_idx numeric. Index of x for random scaling subset inference. Default is 1, the first regressor of x. 
#'    For example, if we want to focus on the 1st and 3rd covariates of x, then set it to be c(1,3).
#' @param level numeric. The confidence level required. Default is 0.95. Can choose 0.90 and 0.80.
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
#' out = sgdi(y~., data=my.dat)

sgdi = function(formula, data, gamma_0=1, alpha=0.667, burn=1, model="lm", inference="rs",
                     bt_start = NULL, qt=0.5,
                     studentize = TRUE, intercept = TRUE, rss_idx = c(1), level=0.95
){
  V_hat_sub = 0
  cl <- match.call()
  
  #----------------------------------------------
  # Linear (Mean) Regression 
  #----------------------------------------------
  if (model=="lm"){
    out = sgdi_lm(formula, data, gamma_0, alpha, burn, inference,
                   bt_start, 
                   studentize, intercept, level
    )
  }
  
  #----------------------------------------------
  # Quantile Regression
  #----------------------------------------------
  if (model=="qr"){
    out = sgdi_qr(formula, data, gamma_0, alpha, burn, inference,
                  bt_start, qt,
                  studentize, intercept, level
    )
  }
  
out$call = cl
return(out)
  
}



