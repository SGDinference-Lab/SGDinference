// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sgd_hqr_cpp
List sgd_hqr_cpp(const arma::mat& x, const arma::colvec& y, const int& burn, const double& gamma_0, const double& alpha, const arma::colvec& bt_start, const double& tau, double h);
RcppExport SEXP _SGDinference_sgd_hqr_cpp(SEXP xSEXP, SEXP ySEXP, SEXP burnSEXP, SEXP gamma_0SEXP, SEXP alphaSEXP, SEXP bt_startSEXP, SEXP tauSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type bt_start(bt_startSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(sgd_hqr_cpp(x, y, burn, gamma_0, alpha, bt_start, tau, h));
    return rcpp_result_gen;
END_RCPP
}
// sgd_qr_cpp
List sgd_qr_cpp(const arma::mat& x, const arma::colvec& y, const int& burn, const double& gamma_0, const double& alpha, const arma::colvec& bt_start, const double& tau);
RcppExport SEXP _SGDinference_sgd_qr_cpp(SEXP xSEXP, SEXP ySEXP, SEXP burnSEXP, SEXP gamma_0SEXP, SEXP alphaSEXP, SEXP bt_startSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type bt_start(bt_startSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(sgd_qr_cpp(x, y, burn, gamma_0, alpha, bt_start, tau));
    return rcpp_result_gen;
END_RCPP
}
// updateGauss
void updateGauss(const arma::mat& Z, const arma::vec& res, arma::vec& der, arma::vec& grad, const double tau, const double h1);
RcppExport SEXP _SGDinference_updateGauss(SEXP ZSEXP, SEXP resSEXP, SEXP derSEXP, SEXP gradSEXP, SEXP tauSEXP, SEXP h1SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type res(resSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type der(derSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grad(gradSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double >::type h1(h1SEXP);
    updateGauss(Z, res, der, grad, tau, h1);
    return R_NilValue;
END_RCPP
}
// sgd_sqr_cpp
List sgd_sqr_cpp(const arma::mat& x, const arma::colvec& y, const int& burn, const double& gamma_0, const double& alpha, const arma::colvec& bt_start, const double& tau, double h);
RcppExport SEXP _SGDinference_sgd_sqr_cpp(SEXP xSEXP, SEXP ySEXP, SEXP burnSEXP, SEXP gamma_0SEXP, SEXP alphaSEXP, SEXP bt_startSEXP, SEXP tauSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type bt_start(bt_startSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(sgd_sqr_cpp(x, y, burn, gamma_0, alpha, bt_start, tau, h));
    return rcpp_result_gen;
END_RCPP
}
// sgdi_lm_cpp
List sgdi_lm_cpp(const arma::mat& x, const arma::colvec& y, const int& burn, const double& gamma_0, const double& alpha, const arma::colvec& bt_start, const std::string inference);
RcppExport SEXP _SGDinference_sgdi_lm_cpp(SEXP xSEXP, SEXP ySEXP, SEXP burnSEXP, SEXP gamma_0SEXP, SEXP alphaSEXP, SEXP bt_startSEXP, SEXP inferenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type bt_start(bt_startSEXP);
    Rcpp::traits::input_parameter< const std::string >::type inference(inferenceSEXP);
    rcpp_result_gen = Rcpp::wrap(sgdi_lm_cpp(x, y, burn, gamma_0, alpha, bt_start, inference));
    return rcpp_result_gen;
END_RCPP
}
// sgdi_qr_cpp
List sgdi_qr_cpp(const arma::mat& x, const arma::colvec& y, const int& burn, const double& gamma_0, const double& alpha, const arma::colvec& bt_start, const std::string inference, const double& tau);
RcppExport SEXP _SGDinference_sgdi_qr_cpp(SEXP xSEXP, SEXP ySEXP, SEXP burnSEXP, SEXP gamma_0SEXP, SEXP alphaSEXP, SEXP bt_startSEXP, SEXP inferenceSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type bt_start(bt_startSEXP);
    Rcpp::traits::input_parameter< const std::string >::type inference(inferenceSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(sgdi_qr_cpp(x, y, burn, gamma_0, alpha, bt_start, inference, tau));
    return rcpp_result_gen;
END_RCPP
}
// sgdi_z_cpp
List sgdi_z_cpp(const arma::mat& x, const arma::colvec& y, const arma::mat& z, const int& burn, const double& gamma_0, const double& alpha, const arma::colvec& bt_start, const std::string inference);
RcppExport SEXP _SGDinference_sgdi_z_cpp(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP burnSEXP, SEXP gamma_0SEXP, SEXP alphaSEXP, SEXP bt_startSEXP, SEXP inferenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const int& >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type bt_start(bt_startSEXP);
    Rcpp::traits::input_parameter< const std::string >::type inference(inferenceSEXP);
    rcpp_result_gen = Rcpp::wrap(sgdi_z_cpp(x, y, z, burn, gamma_0, alpha, bt_start, inference));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SGDinference_sgd_hqr_cpp", (DL_FUNC) &_SGDinference_sgd_hqr_cpp, 8},
    {"_SGDinference_sgd_qr_cpp", (DL_FUNC) &_SGDinference_sgd_qr_cpp, 7},
    {"_SGDinference_updateGauss", (DL_FUNC) &_SGDinference_updateGauss, 6},
    {"_SGDinference_sgd_sqr_cpp", (DL_FUNC) &_SGDinference_sgd_sqr_cpp, 8},
    {"_SGDinference_sgdi_lm_cpp", (DL_FUNC) &_SGDinference_sgdi_lm_cpp, 7},
    {"_SGDinference_sgdi_qr_cpp", (DL_FUNC) &_SGDinference_sgdi_qr_cpp, 8},
    {"_SGDinference_sgdi_z_cpp", (DL_FUNC) &_SGDinference_sgdi_z_cpp, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_SGDinference(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
