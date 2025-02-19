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
NumericVector latlon(NumericVector cd, LogicalVector& value);
RcppExport SEXP _Waypoint_latlon(SEXP cdSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type cd(cdSEXP);
    Rcpp::traits::input_parameter< LogicalVector& >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(latlon(cd, value));
    return rcpp_result_gen;
END_RCPP
}
// printcoords
NumericVector printcoords(NumericVector cd);
RcppExport SEXP _Waypoint_printcoords(SEXP cdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type cd(cdSEXP);
    rcpp_result_gen = Rcpp::wrap(printcoords(cd));
    return rcpp_result_gen;
END_RCPP
}
// validatecoords
NumericVector validatecoords(NumericVector cd);
RcppExport SEXP _Waypoint_validatecoords(SEXP cdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type cd(cdSEXP);
    rcpp_result_gen = Rcpp::wrap(validatecoords(cd));
    return rcpp_result_gen;
END_RCPP
}
// formatcoords
CharacterVector formatcoords(NumericVector nv);
RcppExport SEXP _Waypoint_formatcoords(SEXP nvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nv(nvSEXP);
    rcpp_result_gen = Rcpp::wrap(formatcoords(nv));
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
// printwaypoints
DataFrame printwaypoints(DataFrame df);
RcppExport SEXP _Waypoint_printwaypoints(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(printwaypoints(df));
    return rcpp_result_gen;
END_RCPP
}
// validatewaypoints
DataFrame validatewaypoints(DataFrame df);
RcppExport SEXP _Waypoint_validatewaypoints(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(validatewaypoints(df));
    return rcpp_result_gen;
END_RCPP
}
// as_coord
NumericVector as_coord(DataFrame df, bool latlon);
RcppExport SEXP _Waypoint_as_coord(SEXP dfSEXP, SEXP latlonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< bool >::type latlon(latlonSEXP);
    rcpp_result_gen = Rcpp::wrap(as_coord(df, latlon));
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
    {"_Waypoint_printcoords", (DL_FUNC) &_Waypoint_printcoords, 1},
    {"_Waypoint_validatecoords", (DL_FUNC) &_Waypoint_validatecoords, 1},
    {"_Waypoint_formatcoords", (DL_FUNC) &_Waypoint_formatcoords, 1},
    {"_Waypoint_waypoints", (DL_FUNC) &_Waypoint_waypoints, 2},
    {"_Waypoint_waypoints_replace", (DL_FUNC) &_Waypoint_waypoints_replace, 2},
    {"_Waypoint_printwaypoints", (DL_FUNC) &_Waypoint_printwaypoints, 1},
    {"_Waypoint_validatewaypoints", (DL_FUNC) &_Waypoint_validatewaypoints, 1},
    {"_Waypoint_as_coord", (DL_FUNC) &_Waypoint_as_coord, 2},
    {"_Waypoint_rcpp_hello_world", (DL_FUNC) &_Waypoint_rcpp_hello_world, 0},
    {"_rcpp_module_boot_NumEx", (DL_FUNC) &_rcpp_module_boot_NumEx, 0},
    {"_rcpp_module_boot_yada", (DL_FUNC) &_rcpp_module_boot_yada, 0},
    {"_rcpp_module_boot_stdVector", (DL_FUNC) &_rcpp_module_boot_stdVector, 0},
    {NULL, NULL, 0}
};

void my_package_init(DllInfo *dll);
RcppExport void R_init_Waypoint(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    my_package_init(dll);
}
