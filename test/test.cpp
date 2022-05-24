#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec test(mat &z) {
  int p = 2;
  int q = 5;
  int obs = 4;
  
  Rcout << p << '\n';
  Rcout << q << '\n';
  Rcout << z << '\n';
  
  uvec rnd_sel = randperm(q,p);
  // vec z_rnd = z(obs, rnd_sel);
  
  vec z_rnd = trans(z.row(obs));
  vec z_rnd2 = z_rnd(rnd_sel);
  z_rnd = z_rnd(rnd_sel);
  
  Rcout << rnd_sel << '\n';
  Rcout << z_rnd2 << '\n';

  return z_rnd;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
  
  /*** R
test(matrix(c(1:50),10,5))
*/