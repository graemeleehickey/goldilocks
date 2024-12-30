// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// logrank_instance
std::vector<double> logrank_instance(std::vector<double>& groupa, std::vector<double>& groupb, std::vector<int>& groupacensored, std::vector<int>& groupbcensored, bool onlyz);
RcppExport SEXP _goldilocks_logrank_instance(SEXP groupaSEXP, SEXP groupbSEXP, SEXP groupacensoredSEXP, SEXP groupbcensoredSEXP, SEXP onlyzSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double>& >::type groupa(groupaSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type groupb(groupbSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type groupacensored(groupacensoredSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type groupbcensored(groupbcensoredSEXP);
    Rcpp::traits::input_parameter< bool >::type onlyz(onlyzSEXP);
    rcpp_result_gen = Rcpp::wrap(logrank_instance(groupa, groupb, groupacensored, groupbcensored, onlyz));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_goldilocks_logrank_instance", (DL_FUNC) &_goldilocks_logrank_instance, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_goldilocks(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}