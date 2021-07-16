// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// predictCpp2
NumericVector predictCpp2(NumericVector& xData, NumericVector& yData, NumericVector xPred);
RcppExport SEXP _polijus_predictCpp2(SEXP xDataSEXP, SEXP yDataSEXP, SEXP xPredSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xData(xDataSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type yData(yDataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xPred(xPredSEXP);
    rcpp_result_gen = Rcpp::wrap(predictCpp2(xData, yData, xPred));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_polijus_predictCpp2", (DL_FUNC) &_polijus_predictCpp2, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_polijus(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
