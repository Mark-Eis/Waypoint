# Waypoint (development version)

# Waypoint 1.2.0

* Class and function forward declarations moved to header file CoordBase.h (#113).

* S3 `format()` methods documented more comprehensively (#108).

* Correct error message in `get_coordtype(const int)` (#111).

* Code improved in `format_switch(const T& t)` (#112, #116).

* Remove redundant `Coordbase::get_ff()` (#110).

* Use C++ {fmt} library to ensure formatting and printing of correct widths when names contain extended ASCII codes (#109, #117).

* S3 `format()` and `print()` methods for `"coords"` and `"waypoints"` objects now have a `fmt` argument enabling changing the formatted/printed
  coordinate format (#129, #130, #131). 

# Waypoint 1.1.1

* S3 `format()` method for `"waypoints"` objects `usenames` argument fixed.

* S3 `print()` methods for `"coords"` and `"waypoints"` objects print widths correctly when `max` argument / `getOption("max.print")` is exceeded.

* S3 `validate()` methods for `"coords"` and `"waypoints"` objects now have `force` argument signifying whether to perform full _de novo_ revalidation
  or simply check existing `"valid"`, `"validlat"` and `"validlon"` attributes, essentially to enable the fix to S3 `print()` methods above.

# Waypoint 1.1.0

* Initial CRAN submission.
