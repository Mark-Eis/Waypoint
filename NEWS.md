# Waypoint (development version)

* S3 `format()` methods documented more comprehensively (#108)

* Correct error message in `get_coordtype(const int)` (#111)

* Code improved in `format_switch(const T& t, CoordType ct)` (#112)

* Remove redundant `Coordbase::get_ff()` (#110)

# Waypoint 1.1.1

* S3 `format()` method for `"waypoints"` objects `usenames` argument fixed.

* S3 `print()` methods for `"coords"` and `"waypoints"` objects print widths correctly when `max` argument / `getOption("max.print")` is exceeded.

* Validate S3 methods for `"coords"` and `"waypoints"` objects now have `force` argument signifying whether to perform full _de novo_ revalidation
  or simply check existing `"valid"`, `"validlat"` and `"validlon"` attributes, essentially to enable the fix to S3 `print()` methods above.

# Waypoint 1.1.0

* Initial CRAN submission.
