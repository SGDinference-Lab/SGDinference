#' Census2000
#'
#' The Census2000 dataset 
#' from Acemoglu and Autor (2011) consists of observations on 26,120 nonwhite, female workers. 
#' This small dataset is contructed from "microwage2000_ext.dta" at 
#' \url{https://economics.mit.edu/people/faculty/david-h-autor/data-archive}.
#' Specifically, observations are dropped if hourly wages are missing or 
#' years of education are smaller than 6.    
#' Then, a 5 percent random sample is drawn to make the dataset small.
#' 
#' @references Acemoglu, D. and Autor, D., 2011. 
#' Skills, tasks and technologies: Implications for employment and earnings. 
#' In Handbook of labor economics (Vol. 4, pp. 1043-1171). Elsevier.
#' 
#' @format A data frame with 26,120 rows and 3 variables:
#' \describe{
#' \item{ln_hrwage}{log hourly wages}
#' \item{edyrs}{years of education}
#' \item{exp}{years of potential experience}
#' }
#' 
#' @source The original dataset from Acemoglu and Autor (2011) is available 
#' at \url{https://economics.mit.edu/people/faculty/david-h-autor/data-archive}.
"Census2000"