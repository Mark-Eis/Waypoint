# Changelog

## Waypoint (development version)

- `Coords` as abstract base class (#218).

- Minimised unnecessary dependency on `{fmt}` headers. Dependency on
  `{fmt}` headers is not eliminated completely as it allows
  `stdlenstr(vector<string>&)` to count extended Unicode characters
  (e.g., those with accents) as a single character (#109, \#207).

- Other minor tweaks and tidying of code.

## Waypoint 2.0.0

CRAN release: 2026-07-16

- Extensively revised source code, making use of the newer C++17, C++20
  and C++23 features where possible, for simpler, more understandable
  and more easily maintainable code (#150-#213).

- Static polymorphism and SFINAE used for `FamousFive` stuctures and for
  adding suffixes to formatted coordinates (#195, \#209).

- Fixed edge cases when decimal minutes in `degmin`, or decimal seconds
  in `degminsec`, round up to 60.00 when formatted for printing (#208).

- Shorthand notation for simple, single-type argument concepts instead
  of traditional templated function notation (#175).

- Lambdas for conversion and formatting replaced with specialised
  template function objects managed using generic templated code (#205).

- Type traits and concepts for `NumericVector`, `DataFrame`,
  specialiased vectors and `Coords` class (#186).

- Replace C style `const char*` with `string` (#184).

- Typed vectors to carry coordinate format information, allowing used of
  specialised template functors for conversion and formatting
  (#178-#181, \#203).

- New `Coords` and `Waypoints` classes. Templated `Coords` class has a
  single `NumericVector` representing coordinate values, and `Waypoints`
  has two representing latitude and longitude. `Coords` class implements
  low-level formatting, validation and conversion functions on these
  `NumericVector`s. (#163–#165, \#168, \#171–#174).

- `coordtype_to_int(CoordType)` function adds 1 for consistency with its
  inverse `get_coordtype(int)`, which subtracts 1 (#167).

- Improved and simplified validation algorithms and warnings (#166,
  \#198).

## Waypoint 1.3.0

CRAN release: 2026-06-01

- Improve documentation of `[<-.coords` replacement operator and
  [`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md)
  examples (#160).

- Remove abstruse `convert_switch<>()` function call from
  [`as_coords()`](https://mark-eis.github.io/Waypoint/reference/coords.md)
  and `as_waypoints` (#159).

- Rectify erroneous “Invalid coords!” warning after revalidating valid
  `coords` (#158).

- Simplify `validated()` using bitwise return value and rename as
  `template<NumericVector_or_DataFrame T>`
  `check_logical_attr(T , const char*)` (#157).

- Protect base class functions without public interface (#156).

- Replace `Validator` functor class with lambda in new member function
  `Coordbase::validate0()` (#155).

- Replace `Convert` functor class with lambdas selected using
  `if constexpr … else` statement in new member function
  `Coordbase::convert0()` (#154).

- `if constexpr` statements within templated functions (#153).

- [`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md)
  as pure virtual function in `Coordbase` (#152).

- Replace `Format` functor class with lambdas selected using
  `if constexpr … else` statement in new member function
  `Coordbase::format0()` (#149, \#153).

- Abstract replicated code in `WayPoint::format()` to a single function
  `WayPoint::format2()` (#146).

- `FormatLL<>` functors replaced with lambdas in `Coord::format()` and
  `WayPoint::format()` (#145).

- Replace `static_assert` statements in templates with concepts (#144).

- Remove redundant code from `as_waypoints(DataFrame, int = 1)` (#143).

- Simplify
  `fmt::formatter<CoordType>::format(CoordType, format_context&)`
  (#142).

- Improve `get_coordtype(int i)` (#141, \#142).

## Waypoint 1.2.1

CRAN release: 2025-05-31

- S3 [`print()`](https://rdrr.io/r/base/print.html) methods for
  `"coords"` and `"waypoints"` now employ the null coalescing operator
  `%||%` as intended (#140).

- S3 [`print()`](https://rdrr.io/r/base/print.html) method for
  `"waypoints"` objects now has an explicit `fmt` argument and correct
  formatting of the “Latitude … Longitude” headings when this argument
  is used (#139).

- New S3 extract `` `[`( ``*`<coords>`*`)` and replace
  `` `[<-`( ``*`<coords>`*`)` methods for `"coords"` objects (#135).

- S3 extract `` `[`( ``*`<coords>`*`)` method allows simpler code in
  [`print.coords()`](https://mark-eis.github.io/Waypoint/reference/format.md)
  and
  [`review.coords()`](https://mark-eis.github.io/Waypoint/reference/review.md)
  S3 methods (#136).

- Corrected
  [`as_waypoints()`](https://mark-eis.github.io/Waypoint/reference/waypoints.md)
  and
  [`format()`](https://mark-eis.github.io/Waypoint/reference/format.md)
  documentation (#133, \#134, \#137).

- Note added to documentation for
  [`convert()`](https://mark-eis.github.io/Waypoint/reference/convert.md).

## Waypoint 1.2.0

CRAN release: 2025-05-18

- Class and function forward declarations moved to header file
  CoordBase.h (#113).

- S3
  [`format()`](https://mark-eis.github.io/Waypoint/reference/format.md)
  methods documented more comprehensively (#108).

- Correct error message in `get_coordtype(const int)` (#111).

- Code improved in `format_switch(const T& t)` (#112, \#116).

- Remove redundant `Coordbase::get_ff()` (#110).

- Use C++ {fmt} library to ensure formatting and printing of correct
  widths when names contain extended ASCII codes (#109, \#117).

- S3
  [`format()`](https://mark-eis.github.io/Waypoint/reference/format.md)
  and [`print()`](https://rdrr.io/r/base/print.html) methods for
  `"coords"` and `"waypoints"` objects now have a `fmt` argument
  enabling changing the formatted/printed coordinate format (#129,
  \#130, \#131).

## Waypoint 1.1.1

CRAN release: 2025-04-04

- S3
  [`format()`](https://mark-eis.github.io/Waypoint/reference/format.md)
  method for `"waypoints"` objects `usenames` argument fixed.

- S3 [`print()`](https://rdrr.io/r/base/print.html) methods for
  `"coords"` and `"waypoints"` objects print widths correctly when `max`
  argument / `getOption("max.print")` is exceeded.

- S3
  [`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md)
  methods for `"coords"` and `"waypoints"` objects now have `force`
  argument signifying whether to perform full *de novo* revalidation or
  simply check existing `"valid"`, `"validlat"` and `"validlon"`
  attributes, essentially to enable the fix to S3
  [`print()`](https://rdrr.io/r/base/print.html) methods above.

## Waypoint 1.1.0

CRAN release: 2025-03-19

- Initial CRAN submission.
