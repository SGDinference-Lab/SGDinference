#' @export

plot.sgdi = function(x,path_index=1,..., quote=FALSE) {
#  print(exists("path_index"))
#  if (!exists("path_index")) {
#    path_index=1
#  }
  if ("path_coefficients" %in% ls(x)){
    x_length = nrow(as.matrix(x$path_coefficients))
    plot(x=c(1:x_length), y=x$path_coefficients[,path_index], type="l", xlab="No. of Obs.", ylab="Coefficient Path")
  } else {
    cat("Error: It does not contain 'path_coefficients' to plot. \n")
  }
}