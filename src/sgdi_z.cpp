#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List sgdi_z_cpp(const arma::mat& x, const arma::colvec& y, const arma::mat& z, const int& burn, const double& gamma_0, const double& alpha,
             const arma::colvec& bt_start, const std::string inference){
  int n = y.n_elem;
  double learning_rate_new;
  arma::colvec gradient_bt_new;
  arma::colvec bt_t = bt_start;
  int p = bt_t.n_elem;
  int q = z.n_cols;
  arma::colvec bar_bt_t;
  bar_bt_t.zeros(p);

  arma::mat A_t;
  arma::vec b_t;
  double c_t = 0.0;
  arma::mat V_t;

  A_t.zeros(p,p);
  b_t.zeros(p);
  V_t.zeros(p,p);
  
  double A_t1 = 0.0;
  double b_t1 = 0.0;
  double V_t1 = 0.0;

  if (burn > 0) {
    for(int obs = 1; obs < (burn+1); obs++){
      learning_rate_new = gamma_0 * std::pow(obs, -alpha);
      uvec rnd_sel = randperm(q,p);
      vec z_rnd = trans(z.row(obs-1));
      z_rnd = z_rnd(rnd_sel);
      gradient_bt_new = z_rnd * (x.row(obs-1) * bt_t - y(obs-1));
      bt_t = bt_t - learning_rate_new * gradient_bt_new;
    }
  }

  
/*  rnd_sel = sample.int(q, p, replace=FALSE)
    rnd_z = z[obs, rnd_sel]
    gradient_bt_new = rnd_z %*% (t(x[obs,]) %*% bt_t - y[obs])
*/    
    
    
    
  // for (int obs = burn; obs < (n+1); obs++){
  for (int obs = (burn+1); obs < (n+1); obs++){
    learning_rate_new = gamma_0 * std::pow(obs-burn, -alpha);
    uvec rnd_sel = randperm(q,p);
    vec z_rnd = trans(z.row(obs-1));
    z_rnd = z_rnd(rnd_sel);
    gradient_bt_new = z_rnd * (x.row(obs-1) * bt_t - y(obs-1));
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
                      Named("V_hat") = V_t);
}

