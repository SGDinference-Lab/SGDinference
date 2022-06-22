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

  arma::mat A_t = arma::mat(p,p);
  arma::vec b_t = arma::vec(p);
  double c_t = 0.0;
  arma::mat V_t = arma::mat(p,p);
  
  double A_t1 = 0.0;
  double b_t1 = 0.0;
  double V_t1 = 0.0;

  if (burn > 1) {
    for(int obs = 1; obs < (burn+1); obs++){
      learning_rate_new = gamma_0 * std::pow(obs, -alpha);
      gradient_bt_new = ( trans(x.row(obs-1)) * ( (y(obs-1) < as_scalar(x.row(obs-1) * bt_t)) - tau) );
      bt_t = bt_t - learning_rate_new * gradient_bt_new;
    }
  }

  // for (int obs = burn; obs < (n+1); obs++){
  for (int obs = (burn+1); obs < (n+1); obs++){
    learning_rate_new = gamma_0 * std::pow(obs, -alpha);
    gradient_bt_new = ( trans(x.row(obs-1)) * ( (y(obs-1) < as_scalar(x.row(obs-1) * bt_t)) - tau) );
    bt_t = bt_t - learning_rate_new * gradient_bt_new;
    bar_bt_t = ( bar_bt_t*(obs - burn - 1) + bt_t ) / (obs - burn);
    if ( inference == "rs") {
      A_t = A_t + std::pow(obs - burn, 2.0) * bar_bt_t * trans(bar_bt_t);
      b_t = b_t + std::pow(obs - burn, 2.0) * bar_bt_t;
      c_t = c_t + std::pow(obs - burn, 2.0);
      V_t = ( A_t - b_t * trans(bar_bt_t) - bar_bt_t * trans(b_t) + c_t * bar_bt_t * trans(bar_bt_t) ) / (std::pow(obs - burn, 2.0));
    }
    if ( inference == "rs1") {
      A_t1 = A_t1 + std::pow(obs - burn, 2.0) * bar_bt_t[1] * bar_bt_t[1];
      b_t1 = b_t1 + std::pow(obs - burn, 2.0) * bar_bt_t[1];
      c_t = c_t + std::pow(obs - burn, 2.0);
      V_t1 = ( A_t1 - b_t1 * bar_bt_t[1] - bar_bt_t[1] * b_t1 + c_t * bar_bt_t[1] * bar_bt_t[1] ) / (std::pow(obs - burn, 2.0));
    }
  }


  return List::create(Named("beta_hat") = bar_bt_t,
                      Named("V_hat") = V_t,
                      Named("V_hat1") = V_t1);
}


