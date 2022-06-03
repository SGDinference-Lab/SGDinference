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
#' @param bt_start numeric
#' @param studentize logical. Studentize regressors. Default is TRUE
#' @param intercept logical. Use the intercept term for regressors. Default is TRUE
#' @param qt numeric. Quantile. Default is 0.5. 
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


sgdi = function(x, y, z=NA, gamma_0=1, alpha=0.667, burn=1, model="lm", inference="rs",
                     bt_start = NULL, path_output = NULL, qt=0.5,
                     studentize = TRUE, intercept = TRUE
){
  #----------------------------------------------
  # Linear (Mean) Regression 
  #----------------------------------------------
  if (model=="lm"){
    test = sgdi_lm(x, y, gamma_0, alpha, burn, inference="rs",
                   bt_start, path_output, 
                   studentize, intercept
    )
    beta_hat = test$beta_hat
    V_hat = test$V_hat
  }
  
  
  #----------------------------------------------
  # Quantile Regression
  #----------------------------------------------
  if (model=="qr"){
    out = sgdi_qr(x, y, gamma_0, alpha, burn, inference,
                  bt_start, path_output, qt,
                  studentize, intercept
    )
    beta_hat = out$beta_hat
    V_hat = out$V_hat
  }
  
  if (model=="z"){
    out = sgdi_z(x, y, z, gamma_0, alpha, burn, inference,
                  bt_start, path_output, qt,
                  studentize, intercept
    )
    beta_hat = out$beta_hat
    V_hat = out$V_hat
  }
  
  if ( is.null(path_output)) {
    return(list(beta_hat=beta_hat, V_hat = V_hat))
  } else {
    return(list(beta_hat = beta_hat, V_hat = V_hat, beta_hat_path = beta_hat_path, V_hat_path = V_hat_path))
  }
  
}


# beta_hat_path and V_hat_path are not completed in this code.
