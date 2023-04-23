// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_predvals
arma::mat get_predvals(arma::mat parm_post, arma::mat X_all, int ntotal);
RcppExport SEXP _NBRegAD_get_predvals(SEXP parm_postSEXP, SEXP X_allSEXP, SEXP ntotalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type parm_post(parm_postSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_all(X_allSEXP);
    Rcpp::traits::input_parameter< int >::type ntotal(ntotalSEXP);
    rcpp_result_gen = Rcpp::wrap(get_predvals(parm_post, X_all, ntotal));
    return rcpp_result_gen;
END_RCPP
}
// get_count
arma::mat get_count(arma::vec offset, arma::mat pred_mat, arma::vec overdisp_post);
RcppExport SEXP _NBRegAD_get_count(SEXP offsetSEXP, SEXP pred_matSEXP, SEXP overdisp_postSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pred_mat(pred_matSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type overdisp_post(overdisp_postSEXP);
    rcpp_result_gen = Rcpp::wrap(get_count(offset, pred_mat, overdisp_post));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP _NBRegAD_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NBRegAD_get_predvals", (DL_FUNC) &_NBRegAD_get_predvals, 3},
    {"_NBRegAD_get_count", (DL_FUNC) &_NBRegAD_get_count, 3},
    {"_NBRegAD_rcpp_hello", (DL_FUNC) &_NBRegAD_rcpp_hello, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_NBRegAD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
