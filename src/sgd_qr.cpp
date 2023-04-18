#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List sgd_qr_cpp(const arma::mat& x, const arma::colvec& y, const int& burn, const double& gamma_0, const double& alpha,
             const arma::colvec& bt_start, const double& tau){
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
      gradient_bt_new = ( trans(x.row(obs-1)) * ( (y(obs-1) < as_scalar(x.row(obs-1) * bt_t)) - tau) );
      bt_t = bt_t - learning_rate_new * gradient_bt_new;
    }
  }

  for (int obs = (burn+1); obs < (n+1); obs++){
    learning_rate_new = gamma_0 * std::pow(obs, -alpha);
    gradient_bt_new = ( trans(x.row(obs-1)) * ( (y(obs-1) < as_scalar(x.row(obs-1) * bt_t)) - tau) );
    bt_t = bt_t - learning_rate_new * gradient_bt_new;
    bar_bt_t = ( bar_bt_t*(obs - burn - 1) + bt_t ) / (obs - burn);
  }

  return List::create(Named("beta_hat") = bar_bt_t);
}
