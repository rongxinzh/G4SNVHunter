// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// G4HTranslate
NumericVector G4HTranslate(std::string sequence);
RcppExport SEXP _G4SNVHunter_G4HTranslate(SEXP sequenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type sequence(sequenceSEXP);
    rcpp_result_gen = Rcpp::wrap(G4HTranslate(sequence));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_G4SNVHunter_G4HTranslate", (DL_FUNC) &_G4SNVHunter_G4HTranslate, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_G4SNVHunter(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}