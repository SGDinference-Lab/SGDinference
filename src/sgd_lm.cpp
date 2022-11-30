#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Leading NA
//' 
//' Explanation
//' 
//' @param x numeric. (n x p) matrix of regressors. Should not include 1 (the intercept)
//' @param y numeric
//' @param gamma_0 numeric
//' @param alpha numeric
//' @param burn numeric
//' @param bt_start numeric
//' @export
// [[Rcpp::export]]
List sgd_lm_cpp(const arma::mat& x, const arma::colvec& y, const int& burn, const double& gamma_0, const double& alpha,
             const arma::colvec& bt_start){
  int n = y.n_elem;
  double learning_rate_new;
  arma::colvec gradient_bt_new;
  arma::colvec bt_t = bt_start;
  int p = bt_t.n_elem;
  arma::colvec bar_bt_t;
  bar_bt_t.zeros(p);

  if (burn > 1) {
    for(int obs = 1; obs < (burn+1); obs++){
      learning_rate_new = gamma_0 * std::pow(obs, -alpha);
      gradient_bt_new = trans(x.row(obs-1)) * (x.row(obs-1) * bt_t - y(obs-1));
      bt_t = bt_t - learning_rate_new * gradient_bt_new;
    }
  }

  for (int obs = (burn+1); obs < (n+1); obs++){
    learning_rate_new = gamma_0 * std::pow(obs, -alpha);
    gradient_bt_new = trans(x.row(obs-1)) * (x.row(obs-1) * bt_t - y(obs-1));
    bt_t = bt_t - learning_rate_new * gradient_bt_new;
    bar_bt_t = ( bar_bt_t*(obs - burn - 1) + bt_t ) / (obs - burn);
  }
  

  //-------------------------------------------

  return List::create(Named("beta_hat") = bar_bt_t);
}

