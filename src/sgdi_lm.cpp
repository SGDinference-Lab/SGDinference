#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List sgdi_lm_cpp(const arma::mat& x, 
                 const arma::colvec& y, 
                 const int& burn, 
                 const double& gamma_0, 
                 const double& alpha,
                 const arma::colvec& bt_start, 
                 const std::string inference,
                 const arma::uvec& rss_idx,
                 const arma::rowvec& x_mean,
                 const arma::rowvec& x_sd,
                 const bool& path,
                 const arma::colvec& path_index) {
  int n = y.n_elem;
  double learning_rate_new;
  arma::colvec gradient_bt_new;
  arma::colvec bt_t = bt_start;
  int p = bt_t.n_elem;
  int n_s;
  arma::colvec bar_bt_t = bt_start;
  
  // Initialize the path output variables when {path} is true
  int p_path = path_index.n_elem;
  arma::mat bar_bt_path = arma::mat(n,p_path);
  for (int i_p_path = 0; i_p_path < p_path; i_p_path++) {
    bar_bt_path(0, i_p_path) = bt_start(path_index(i_p_path)-1);
  }
  
  arma::mat A_t = arma::mat(p,p);
  arma::vec b_t = arma::vec(p);
  double c_t = 0.0;
  arma::mat V_t = arma::mat(p,p);
  
  // Random scaling subset inference: arma::matrix initialization
  int p_s = rss_idx.n_elem;
  arma::mat A_ts = arma::mat(p_s,p_s);
  arma::vec b_ts = arma::vec(p_s);
  arma::mat V_ts = arma::mat(p_s,p_s);
  
  // Random scaling diagonal inference: arma::matrix initialization
  arma::vec A_td = arma::vec(p);
  arma::vec b_td = arma::vec(p);
  arma::vec V_td = arma::vec(p);
  
  if (x_sd(0) < 0) {
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
      n_s = obs - burn;
      bar_bt_t = ( bar_bt_t*(n_s - 1) + bt_t ) / (n_s);
      
      // Save the bar_bt_t for the path outcome
      if (path) {
        for (int i_p_path = 0; i_p_path < p_path; i_p_path++) {
          bar_bt_path(obs-1, i_p_path) = bar_bt_t(path_index(i_p_path)-1);
        }
      }
      
      if ( inference == "rs") {
        A_t = A_t + std::pow(n_s, 2.0) * bar_bt_t * trans(bar_bt_t);
        b_t = b_t + std::pow(n_s, 2.0) * bar_bt_t;
        c_t = c_t + std::pow(n_s, 2.0);
      }
      if ( inference == "rss") {
        A_ts = A_ts + std::pow(n_s, 2.0) * bar_bt_t(rss_idx) * trans(bar_bt_t(rss_idx));
        b_ts = b_ts + std::pow(n_s, 2.0) * bar_bt_t(rss_idx);
        c_t = c_t + std::pow(n_s, 2.0);
      }
      if ( inference == "rsd") {
        A_td = A_td + std::pow(n_s, 2.0) * (bar_bt_t % bar_bt_t);
        b_td = b_td + std::pow(n_s, 2.0) * bar_bt_t;
        c_t = c_t + std::pow(n_s, 2.0);      V_td = ( A_td - 2*( b_td % bar_bt_t ) + c_t * ( bar_bt_t % bar_bt_t ) ) / (std::pow(n_s, 2.0));
      }
    }
  } else {
    if (burn > 1) {
      for(int obs = 1; obs < (burn+1); obs++){
        learning_rate_new = gamma_0 * std::pow(obs, -alpha);
        gradient_bt_new = trans((x.row(obs-1)-x_mean)/x_sd) * ((x.row(obs-1)-x_mean)/x_sd * bt_t - y(obs-1));
        bt_t = bt_t - learning_rate_new * gradient_bt_new;
      }
    }
  
    for (int obs = (burn+1); obs < (n+1); obs++){
      learning_rate_new = gamma_0 * std::pow(obs, -alpha);
      gradient_bt_new = trans((x.row(obs-1)-x_mean)/x_sd) * ((x.row(obs-1)-x_mean)/x_sd * bt_t - y(obs-1));
      bt_t = bt_t - learning_rate_new * gradient_bt_new;
      n_s = obs - burn;
      bar_bt_t = ( bar_bt_t*(n_s - 1) + bt_t ) / (n_s);
      
      // Save the bar_bt_t for the path outcome
      if (path) {
        for (int i_p_path = 0; i_p_path < p_path; i_p_path++) {
          bar_bt_path(obs-1, i_p_path) = bar_bt_t(path_index(i_p_path)-1);
        }
      }
      
      if ( inference == "rs") {
        A_t = A_t + std::pow(n_s, 2.0) * bar_bt_t * trans(bar_bt_t);
        b_t = b_t + std::pow(n_s, 2.0) * bar_bt_t;
        c_t = c_t + std::pow(n_s, 2.0);
      }
      if ( inference == "rss") {
        A_ts = A_ts + std::pow(n_s, 2.0) * bar_bt_t(rss_idx) * trans(bar_bt_t(rss_idx));
        b_ts = b_ts + std::pow(n_s, 2.0) * bar_bt_t(rss_idx);
        c_t = c_t + std::pow(n_s, 2.0);
      }
      if ( inference == "rsd") {
        A_td = A_td + std::pow(n_s, 2.0) * (bar_bt_t % bar_bt_t);
        b_td = b_td + std::pow(n_s, 2.0) * bar_bt_t;
        c_t = c_t + std::pow(n_s, 2.0);      V_td = ( A_td - 2*( b_td % bar_bt_t ) + c_t * ( bar_bt_t % bar_bt_t ) ) / (std::pow(n_s, 2.0));
      }
    }
  }
  
  n_s = n - burn;
  if ( inference == "rs") {
    V_t = ( A_t - b_t * trans(bar_bt_t) - bar_bt_t * trans(b_t) + c_t * bar_bt_t * trans(bar_bt_t) ) / (std::pow(n_s, 2.0));
  }
  if ( inference == "rss") {
    V_ts = ( A_ts - b_ts * trans(bar_bt_t(rss_idx)) - bar_bt_t(rss_idx) * trans(b_ts) + c_t * bar_bt_t(rss_idx) * trans(bar_bt_t(rss_idx)) ) / (std::pow(n_s, 2.0));
  }
  if ( inference == "rsd") {
    V_td = ( A_td - 2*( b_td % bar_bt_t ) + c_t * ( bar_bt_t % bar_bt_t ) ) / (std::pow(n_s, 2.0));
  }
  
  return List::create(Named("beta_hat") = bar_bt_t,
                      Named("V_hat") = V_t,
                      Named("V_hat_sub") = V_ts,
                      Named("V_hat_diag") = V_td,
                      Named("beta_hat_path") = bar_bt_path);
}