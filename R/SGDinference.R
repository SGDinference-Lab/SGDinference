#' SGDinference
#' 
#' The SGDinference package provides estimation and inference methods for large-scale mean and quantile regression models via stochastic (sub-)gradient descent (S-subGD) algorithms. 
#' The inference procedure handles cross-sectional data sequentially: 
#' (i) updating the parameter estimate with each incoming "new observation", 
#' (ii) aggregating it as a Polyak-Ruppert average, and 
#' (iii) computing an asymptotically pivotal statistic for inference through random scaling.
#' 
#' @docType package
#' @author Sokbae Lee, Yuan Liao, Myung Hwan Seo, Youngki Shin
#' @import Rcpp stats
#' @importFrom Rcpp sourceCpp
#' @useDynLib SGDinference
#' @name SGDinference
NULL
