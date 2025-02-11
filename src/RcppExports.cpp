// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// coords
NumericVector coords(NumericVector nv, const int fmt);
RcppExport SEXP _Waypoint_coords(SEXP nvSEXP, SEXP fmtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nv(nvSEXP);
    Rcpp::traits::input_parameter< const int >::type fmt(fmtSEXP);
    rcpp_result_gen = Rcpp::wrap(coords(nv, fmt));
    return rcpp_result_gen;
END_RCPP
}
// coords_replace
NumericVector coords_replace(NumericVector nv, int value);
RcppExport SEXP _Waypoint_coords_replace(SEXP nvSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nv(nvSEXP);
    Rcpp::traits::input_parameter< int >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(coords_replace(nv, value));
    return rcpp_result_gen;
END_RCPP
}
// latlon
NumericVector latlon(NumericVector nv, LogicalVector& value);
RcppExport SEXP _Waypoint_latlon(SEXP nvSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nv(nvSEXP);
    Rcpp::traits::input_parameter< LogicalVector& >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(latlon(nv, value));
    return rcpp_result_gen;
END_RCPP
}
// printcoord
NumericVector printcoord(NumericVector nv);
RcppExport SEXP _Waypoint_printcoord(SEXP nvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nv(nvSEXP);
    rcpp_result_gen = Rcpp::wrap(printcoord(nv));
    return rcpp_result_gen;
END_RCPP
}
// validatecoord
NumericVector validatecoord(NumericVector nv);
RcppExport SEXP _Waypoint_validatecoord(SEXP nvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nv(nvSEXP);
    rcpp_result_gen = Rcpp::wrap(validatecoord(nv));
    return rcpp_result_gen;
END_RCPP
}
// formatcoord
CharacterVector formatcoord(NumericVector nv);
RcppExport SEXP _Waypoint_formatcoord(SEXP nvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nv(nvSEXP);
    rcpp_result_gen = Rcpp::wrap(formatcoord(nv));
    return rcpp_result_gen;
END_RCPP
}
// waypoints
DataFrame waypoints(DataFrame df, int fmt);
RcppExport SEXP _Waypoint_waypoints(SEXP dfSEXP, SEXP fmtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type fmt(fmtSEXP);
    rcpp_result_gen = Rcpp::wrap(waypoints(df, fmt));
    return rcpp_result_gen;
END_RCPP
}
// waypoints_replace
DataFrame waypoints_replace(DataFrame df, int value);
RcppExport SEXP _Waypoint_waypoints_replace(SEXP dfSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(waypoints_replace(df, value));
    return rcpp_result_gen;
END_RCPP
}
// printwaypoint
DataFrame printwaypoint(DataFrame df);
RcppExport SEXP _Waypoint_printwaypoint(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(printwaypoint(df));
    return rcpp_result_gen;
END_RCPP
}
// validatewaypoint
DataFrame validatewaypoint(DataFrame df);
RcppExport SEXP _Waypoint_validatewaypoint(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(validatewaypoint(df));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _Waypoint_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_NumEx();
RcppExport SEXP _rcpp_module_boot_yada();
RcppExport SEXP _rcpp_module_boot_stdVector();

static const R_CallMethodDef CallEntries[] = {
    {"_Waypoint_coords", (DL_FUNC) &_Waypoint_coords, 2},
    {"_Waypoint_coords_replace", (DL_FUNC) &_Waypoint_coords_replace, 2},
    {"_Waypoint_latlon", (DL_FUNC) &_Waypoint_latlon, 2},
    {"_Waypoint_printcoord", (DL_FUNC) &_Waypoint_printcoord, 1},
    {"_Waypoint_validatecoord", (DL_FUNC) &_Waypoint_validatecoord, 1},
    {"_Waypoint_formatcoord", (DL_FUNC) &_Waypoint_formatcoord, 1},
    {"_Waypoint_waypoints", (DL_FUNC) &_Waypoint_waypoints, 2},
    {"_Waypoint_waypoints_replace", (DL_FUNC) &_Waypoint_waypoints_replace, 2},
    {"_Waypoint_printwaypoint", (DL_FUNC) &_Waypoint_printwaypoint, 1},
    {"_Waypoint_validatewaypoint", (DL_FUNC) &_Waypoint_validatewaypoint, 1},
    {"_Waypoint_rcpp_hello_world", (DL_FUNC) &_Waypoint_rcpp_hello_world, 0},
    {"_rcpp_module_boot_NumEx", (DL_FUNC) &_rcpp_module_boot_NumEx, 0},
    {"_rcpp_module_boot_yada", (DL_FUNC) &_rcpp_module_boot_yada, 0},
    {"_rcpp_module_boot_stdVector", (DL_FUNC) &_rcpp_module_boot_stdVector, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_Waypoint(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
