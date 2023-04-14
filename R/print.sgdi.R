#' @export

print.sgdi = function(x, ..., quote=FALSE) {
  cat("Call: \n")
  print(x$call, quote=quote)
  cat("\n")
  cat("Estimation Results: \n")
  if (is.null(x$ci.lower)){
    inference.out = data.frame(Coefficient=x$coefficient)    
    rownames(inference.out) = c("Intercept", attr(x$terms,"term.labels"))
    print(inference.out, quote=quote) 
  } else {
    inference.out = data.frame(Coefficient=x$coefficient, CI.Lower=x$ci.lower, CI.Upper = x$ci.upper)    
    rownames(inference.out) = c("Intercept", attr(x$terms,"term.labels"))
    if (x$inference == "rsd" ){
      print(inference.out[-1,], quote=quote)
    } else if (x$inference == "rss"){
      print(inference.out[x$rss_idx_r,], quote=quote)
    } else {
      print(inference.out, quote=quote) 
    }
  }
}