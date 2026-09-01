# Waypoint 2.0.1

* Initialise lambda captures with expressions, and use `mutable` specifier (#221).

* Make use of (C++20 / C++23) Ranges library (#220).

* `Coords` made an abstract base class (#218).

* Simplify `Coords` class constructor (#215).

* Improved messages in `convertcoords()` and `convertwaypoints()` when `fmt` unchanged (#214).

* Minimised dependency on {fmt} headers in favour of conventional `<iostream>` and `<sstream>`
  headers. Dependency not eliminated entirely; {fmt} headers remain essential to enable
  `stdlenstr(vector<string>&)` to count extended Unicode characters (e.g., those with accents)
  as single, rather than as two characters (#109, #216, #217).

* Other minor tweaks and tidying of code.



# Waypoint 2.0.0

* Extensively revised source code, making use of the newer C++17, C++20 and C++23 features where
  possible, for simpler, more understandable and more easily maintainable code (#150-#213).
  
* Static polymorphism and SFINAE used for `FamousFive` stuctures and for adding suffixes to formatted
  coordinates (#195, #209).

* Fixed edge cases when decimal minutes in `degmin`, or decimal seconds in `degminsec`, round up
  to 60.00 when formatted for printing (#208).

* Shorthand notation for simple, single-type argument concepts instead of traditional templated
  function notation (#175).
  
* Lambdas for conversion and formatting replaced with specialised template function objects managed
  using generic templated code (#205).

* Type traits and concepts for `NumericVector`, `DataFrame`, specialiased vectors and `Coords`
  class (#186).

* Replace C style `const char*` with `string` (#184).

* Typed vectors to carry coordinate format information, allowing used of specialised template
  functors for conversion and formatting (#178-#181, #203).

* New `Coords` and `Waypoints` classes. Templated `Coords` class has a single `NumericVector`
  representing coordinate values, and `Waypoints` has two representing latitude and longitude.
  `Coords` class implements low-level formatting, validation and conversion functions on these
  `NumericVector`s. (#163–#165, #168, #171–#174).

* `coordtype_to_int(CoordType)` function adds 1 for consistency with its inverse `get_coordtype(int)`,
  which subtracts 1 (#167).
   
* Improved and simplified validation algorithms and warnings (#166, #198).

# Waypoint 1.3.0

* Improve documentation of `[<-.coords` replacement operator and `validate()` examples (#160).

* Remove abstruse `convert_switch<>()` function call from `as_coords()` and `as_waypoints` (#159).

* Rectify  erroneous "Invalid coords!" warning after revalidating valid `coords` (#158).

* Simplify `validated()` using bitwise return value and rename as `template<NumericVector_or_DataFrame T>`
  `check_logical_attr(T , const char*)` (#157).

* Protect base class functions without public interface (#156).

* Replace `Validator` functor class with lambda in new member function `Coordbase::validate0()` (#155).

* Replace `Convert` functor class with lambdas selected using `if constexpr … else` statement in new
  member function `Coordbase::convert0()` (#154).

* `if constexpr` statements within templated functions (#153).

* `validate()` as pure virtual function in `Coordbase` (#152).

* Replace `Format` functor class with lambdas selected using `if constexpr … else` statement in new
  member function `Coordbase::format0()` (#149, #153).

* Abstract replicated code in `WayPoint::format()` to a single function `WayPoint::format2()` (#146).

* `FormatLL<>` functors replaced with lambdas in `Coord::format()` and `WayPoint::format()` (#145).

* Replace `static_assert` statements in templates with concepts (#144).

* Remove redundant code from `as_waypoints(DataFrame, int = 1)` (#143).

* Simplify `fmt::formatter<CoordType>::format(CoordType, format_context&)` (#142).

* Improve `get_coordtype(int i)` (#141, #142).

# Waypoint 1.2.1

* S3 `print()` methods for `"coords"` and `"waypoints"` now employ the null coalescing operator
  `%||%` as intended (#140). 

* S3 `print()` method for `"waypoints"` objects now has an explicit `fmt` argument and correct
  formatting of the "Latitude ... Longitude" headings when this argument is used (#139).

* New S3 extract <code>&#96;[&#96;(<i>&lt;coords&gt;</i>)</code> and replace
  <code>&#96;[<-&#96;(<i>&lt;coords&gt;</i>)</code> methods for `"coords"` objects (#135).

* S3 extract <code>&#96;[&#96;(<i>&lt;coords&gt;</i>)</code> method allows simpler code in
  `print.coords()` and `review.coords()` S3 methods (#136).

* Corrected `as_waypoints()` and `format()` documentation (#133, #134, #137).

* Note added to documentation for `convert()`.

# Waypoint 1.2.0

* Class and function forward declarations moved to header file CoordBase.h (#113).

* S3 `format()` methods documented more comprehensively (#108).

* Correct error message in `get_coordtype(const int)` (#111).

* Code improved in `format_switch(const T& t)` (#112, #116).

* Remove redundant `Coordbase::get_ff()` (#110).

* Use C++ {fmt} library to ensure formatting and printing of correct widths when names contain extended ASCII
  codes (#109, #117).

* S3 `format()` and `print()` methods for `"coords"` and `"waypoints"` objects now have a `fmt` argument enabling
  changing the formatted/printed coordinate format (#129, #130, #131). 

# Waypoint 1.1.1

* S3 `format()` method for `"waypoints"` objects `usenames` argument fixed.

* S3 `print()` methods for `"coords"` and `"waypoints"` objects print widths correctly when `max` argument / `getOption("max.print")` is exceeded.

* S3 `validate()` methods for `"coords"` and `"waypoints"` objects now have `force` argument signifying whether
  to perform full _de novo_ revalidation or simply check existing `"valid"`, `"validlat"` and `"validlon"`
  attributes, essentially to enable the fix to S3 `print()` methods above.

# Waypoint 1.1.0

* Initial CRAN submission.
