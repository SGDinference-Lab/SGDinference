#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List sgd_hqr_cpp(const arma::mat& x, const arma::colvec& y, const int& burn, const double& gamma_0, const double& alpha,
             const arma::colvec& bt_start, const double& tau, double h=0.05){
  int n = y.n_elem;
  double learning_rate_new;
  arma::colvec gradient_bt_new;
  arma::colvec bt_t = bt_start;
  int p = bt_t.n_elem;
  arma::colvec bar_bt_t;
  bar_bt_t.zeros(p);
  
  if (h <= 0.05) {
    h = std::max(std::pow((std::log(n) + p) / n, 0.4), 0.05);
  }
  const double h1 = 1.0 / h;

  if (burn > 1) {
    for(int obs = 1; obs < burn; obs++){
      learning_rate_new = gamma_0 * std::pow(obs, -alpha);
      double arg =  (-1.0) * ( y(obs-1) - as_scalar(x.row(obs-1)*bt_t) ) * h1;
      gradient_bt_new =  trans(x.row(obs-1)) * ( ( arma::normcdf(arg) - tau ) + arg * arma::normpdf(arg));
      bt_t = bt_t - learning_rate_new * gradient_bt_new;
    }
  }



  // for (int obs = burn; obs < (n+1); obs++){
  for (int obs = burn; obs < (n+1); obs++){
    learning_rate_new = gamma_0 * std::pow(obs, -alpha);
    double arg =  (-1.0) * ( y(obs-1) - as_scalar(x.row(obs-1)*bt_t) ) * h1;
    gradient_bt_new =  trans(x.row(obs-1)) * ( ( arma::normcdf(arg) - tau ) + arg * arma::normpdf(arg));
    bt_t = bt_t - learning_rate_new * gradient_bt_new;
    bar_bt_t = ( bar_bt_t*(obs - burn) + bt_t ) / (obs - burn + 1);
  }


  return List::create(Named("beta_hat") = bar_bt_t);
}


