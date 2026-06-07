# Waypoint (development version)

* `coordtype_to_int(CoordType)` adds 1 for consistency with its inverse function `get_coordtype(int)`,
   which subtracts 1 (#167).

* `template<NumericVector_or_DataFrame T, Coords_or_Waypoints U> const T validate(const T t)`, absorbed
  into `template<NumericVector_or_DataFrame T, Coords_or_Waypoints U> bool revalidate(const T)` (#166).

* `prefixvecstr(vector<string>&, const vector<T>&)` simply overloaded rather than templated (#162).

* Introduce enum class `CoordType` type traits (#161).

* `FamousFive` classes now combine generic and OO techniques in abstract non-template base class with pure
   virtual functions inherited in with three templated derived classes, and instantiated in each `Cordlet`
   class object (#151).
   
*  Improved and simplified validation warnings.

*  New `Coords` and `Waypoints` classes each inheriting from abstract base class `CrdWptBase`, which implements
   member functions common to both derived classes or as pure virtual functions where the two derived classes
   differ. Class `Coords` has a single `NumericVector` representing coordinate values, and `Waypoints` has two
   representing latitude and longitude. `Coordlet` class implements low-level formatting, validation and
   conversion functions on these `NumericVector`s.

*  Extensively revised source code, making use of the newer C++17, C++20 and C++23 feature where possible for
   simpler, more understandable and maintainable code .

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
