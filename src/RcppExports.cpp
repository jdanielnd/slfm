// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// gibbs
List gibbs(NumericMatrix x, int ite, double a, double b, double gamma_a, double gamma_b, double omega_0, double omega_1, bool degenerate);
RcppExport SEXP _slfm_gibbs(SEXP xSEXP, SEXP iteSEXP, SEXP aSEXP, SEXP bSEXP, SEXP gamma_aSEXP, SEXP gamma_bSEXP, SEXP omega_0SEXP, SEXP omega_1SEXP, SEXP degenerateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ite(iteSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_a(gamma_aSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_b(gamma_bSEXP);
    Rcpp::traits::input_parameter< double >::type omega_0(omega_0SEXP);
    Rcpp::traits::input_parameter< double >::type omega_1(omega_1SEXP);
    Rcpp::traits::input_parameter< bool >::type degenerate(degenerateSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs(x, ite, a, b, gamma_a, gamma_b, omega_0, omega_1, degenerate));
    return rcpp_result_gen;
END_RCPP
}
// slfm_MDN
List slfm_MDN(NumericMatrix x, double a, double b, double gamma_a, double gamma_b, double omega_1, int burnin, int lag, int npost);
RcppExport SEXP _slfm_slfm_MDN(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP gamma_aSEXP, SEXP gamma_bSEXP, SEXP omega_1SEXP, SEXP burninSEXP, SEXP lagSEXP, SEXP npostSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_a(gamma_aSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_b(gamma_bSEXP);
    Rcpp::traits::input_parameter< double >::type omega_1(omega_1SEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< int >::type npost(npostSEXP);
    rcpp_result_gen = Rcpp::wrap(slfm_MDN(x, a, b, gamma_a, gamma_b, omega_1, burnin, lag, npost));
    return rcpp_result_gen;
END_RCPP
}
// slfm_MNN
List slfm_MNN(NumericMatrix x, double a, double b, double gamma_a, double gamma_b, double omega_0, double omega_1, int burnin, int lag, int npost);
RcppExport SEXP _slfm_slfm_MNN(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP gamma_aSEXP, SEXP gamma_bSEXP, SEXP omega_0SEXP, SEXP omega_1SEXP, SEXP burninSEXP, SEXP lagSEXP, SEXP npostSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_a(gamma_aSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_b(gamma_bSEXP);
    Rcpp::traits::input_parameter< double >::type omega_0(omega_0SEXP);
    Rcpp::traits::input_parameter< double >::type omega_1(omega_1SEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< int >::type npost(npostSEXP);
    rcpp_result_gen = Rcpp::wrap(slfm_MNN(x, a, b, gamma_a, gamma_b, omega_0, omega_1, burnin, lag, npost));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_slfm_gibbs", (DL_FUNC) &_slfm_gibbs, 9},
    {"_slfm_slfm_MDN", (DL_FUNC) &_slfm_slfm_MDN, 9},
    {"_slfm_slfm_MNN", (DL_FUNC) &_slfm_slfm_MNN, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_slfm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
