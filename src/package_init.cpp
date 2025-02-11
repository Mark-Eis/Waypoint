#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::init]]
void my_package_init(DllInfo *dll) {
	Rcout << "\n*** Initialising Waypoint Package (MCE 2025) ***\n\n";
}