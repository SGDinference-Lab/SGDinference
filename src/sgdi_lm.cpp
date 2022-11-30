#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List sgdi_lm_cpp(const arma::mat& x, const arma::colvec& y, const int& burn, const double& gamma_0, const double& alpha,
             const arma::colvec& bt_start, const std::string inference,
             const arma::uvec& rss_idx){
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
  
  // Random scaling subset inference: arma::matrix initialization
  int p_s = rss_idx.n_elem;
  arma::mat A_ts = arma::mat(p_s,p_s);
  arma::vec b_ts = arma::vec(p_s);
  arma::mat V_ts = arma::mat(p_s,p_s);
  
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
    if ( inference == "rs") {
      A_t = A_t + std::pow(obs - burn, 2.0) * bar_bt_t * trans(bar_bt_t);
      b_t = b_t + std::pow(obs - burn, 2.0) * bar_bt_t;
      c_t = c_t + std::pow(obs - burn, 2.0);
      V_t = ( A_t - b_t * trans(bar_bt_t) - bar_bt_t * trans(b_t) + c_t * bar_bt_t * trans(bar_bt_t) ) / (std::pow(obs - burn, 2.0));
    }
    if ( inference == "rss") {
      A_ts = A_ts + std::pow(obs - burn, 2.0) * bar_bt_t(rss_idx) * trans(bar_bt_t(rss_idx));
      b_ts = b_ts + std::pow(obs - burn, 2.0) * bar_bt_t(rss_idx);
      c_t = c_t + std::pow(obs - burn, 2.0);
      V_ts = ( A_ts - b_ts * trans(bar_bt_t(rss_idx)) - bar_bt_t(rss_idx) * trans(b_ts) + c_t * bar_bt_t(rss_idx) * trans(bar_bt_t(rss_idx)) ) / (std::pow(obs - burn, 2.0));
    }
  }
  
  return List::create(Named("beta_hat") = bar_bt_t,
                      Named("V_hat") = V_t,
                      Named("V_hat_sub") = V_ts);
}