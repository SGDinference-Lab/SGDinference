#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List sgdi_qr_cpp(const arma::mat& x, const arma::colvec& y, const int& burn, const double& gamma_0, const double& alpha,
             const arma::colvec& bt_start, const std::string inference, const double& tau){
  int n = y.n_elem;
  double learning_rate_new;
  arma::colvec gradient_bt_new;
  arma::colvec bt_t = bt_start;
  int p = bt_t.n_elem;
  arma::colvec bar_bt_t;
  bar_bt_t.zeros(p);

  arma::mat A_t;
  arma::vec b_t;
  double c_t = 0.0;
  arma::mat V_t;

  A_t.zeros(p,p);
  b_t.zeros(p);
  V_t.zeros(p,p);

  if (burn > 1) {
    for(int obs = 1; obs < burn; obs++){
      learning_rate_new = gamma_0 * std::pow(obs, -alpha);
      gradient_bt_new = ( trans(x.row(obs-1)) * ( (y(obs-1) < as_scalar(x.row(obs-1) * bt_t)) - tau) );
      bt_t = bt_t - learning_rate_new * gradient_bt_new;
    }
  }



  // for (int obs = burn; obs < (n+1); obs++){
  for (int obs = burn; obs < (n+1); obs++){
    learning_rate_new = gamma_0 * std::pow(obs, -alpha);
    gradient_bt_new = ( trans(x.row(obs-1)) * ( (y(obs-1) < as_scalar(x.row(obs-1) * bt_t)) - tau) );
    bt_t = bt_t - learning_rate_new * gradient_bt_new;
    bar_bt_t = ( bar_bt_t*(obs - burn) + bt_t ) / (obs - burn + 1);
    if ( inference == "rs") {
      A_t = A_t + std::pow(obs, 2.0) * bar_bt_t * trans(bar_bt_t);
      b_t = b_t + std::pow(obs, 2.0) * bar_bt_t;
      c_t = c_t + std::pow(obs, 2.0);
      V_t = ( A_t - b_t * trans(bar_bt_t) - bar_bt_t * trans(b_t) + c_t * bar_bt_t * trans(bar_bt_t) ) / (std::pow(obs, 2.0));
    }
  }


  return List::create(Named("beta_hat") = bar_bt_t,
                      Named("V_hat") = V_t);
}


