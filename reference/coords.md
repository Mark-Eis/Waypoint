# Geographic Coordinate Class

`as_coords()` creates an object of class `"coords"`, a robust
representation of a series of geographic or GPS coordinate values.

`latlon()<-` adds information to objects of class `"coords"` specifying
whether individual coordinate values represent latitude or longitude.

## Usage

``` r
as_coords(object, ...)

# Default S3 method
as_coords(object, ..., fmt = 1L)

latlon(cd) <- value

# S3 method for class 'waypoints'
as_coords(object, which, ...)
```

## Arguments

- object:

  a `numeric` vector of coordinate values, optionally named, or an
  object of class `"waypoints"`.

- ...:

  further arguments passed to or from other methods.

- fmt:

  `integer`, `1L`, `2L` or `3L`, specifying the required coordinate
  format.

- cd:

  object of class `"coords"` created by function `as_coords()`.

- value:

  a `logical` vector of length `1` or `length(x)`.

- which:

  `logical`, indicating whether the `as_coords()` method for class
  `"waypoints"` extracts the latitude component of argument `object` (if
  `TRUE`), or the longitude (if `FALSE`).

## Value

`as_cords()` returns an object of class `"coords"`, comprising a
`numeric` vector argument with additional attributes: –

- `"class"`:

  the `character` string "coords".

- `"fmt"`:

  an `integer` representing the coordinate format.

- `"valid"`:

  a `logical` vector indicating whether individual coordinate values are
  valid geographic locations.

The `as_cords()` method for class `"coords"` returns its `numeric`
vector argument `object` modified in place, whereas the method for class
'waypoints' returns a new `numeric` vector.

`latlon()<-` returns its `"coords"` argument `cd` with a `logical`
vector attribute `"latlon"` added or updated to reflect argument
`value`.

## Details

Individual values provided in a `numeric` vector argument `object`
should have a decimal point after the number of whole degrees in the
case of *decimal degrees*, after the number of whole minutes in the case
of *degrees and minutes*, and after the number of whole seconds in the
case of *degrees, minutes and seconds*.

The `fmt` argument should be `1L` to represent decimal degrees, `2L` for
degrees and minutes, and `3L` for degrees, minutes and seconds and is
used to provide the format of values in the `numeric` vector argument
`object` to be converted to class `"coords"`.

The values of a newly created `"coords"` object are checked to ensure
they are valid geographic locations as described under
[`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md).

Individual coordinate values in a `Coords` object may be specified as
representing latitude or longitude using `latlon()<-`. The `value`
argument may either be a single value, `TRUE` signifying that all values
are latitude, `FALSE` signifying that all values are longitude, or a
`logical` vector of the same length as as the `Coords` object signifying
whether individual values are latitude or longitude.

## See also

[`attr()`](https://rdrr.io/r/base/attr.html),
[`attributes`](https://rdrr.io/r/base/attributes.html),
[`format`](https://mark-eis.github.io/Waypoint/reference/format.md)`(<coords>)`,
[`print`](https://mark-eis.github.io/Waypoint/reference/format.md)`(<coords>)`,
[`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md)
and
[`[<-`](https://mark-eis.github.io/Waypoint/reference/Extract.md)`(<coords>)`.

Other coordsandway:
[`convert()`](https://mark-eis.github.io/Waypoint/reference/convert.md),
[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.md)

## Examples

``` r
## Numeric vector representing degrees and minutes, with
## the decimal point after the number of whole minutes
dm <- c(5130.4659, 4932.7726, 4806.4339, 3853.3696, 0.0000, -3706.7044, -5306.2869, -2514.4093,
        -007.6754, 1823.9137, -12246.7203, -7702.1145, 0.0000, -1217.3178, 7331.0370, -5731.1536)

## Create an unnamed "coords" object in degrees and minutes (fmt = 2)
## (Latitude and longitude unspecified)
as_coords(dm, fmt = 2)
#>  51°30.4659′ (N/E)
#>  49°32.7726′ (N/E)
#>  48°06.4339′ (N/E)
#>  38°53.3696′ (N/E)
#>   0°00.0000′ (N/E)
#>  37°06.7044′ (S/W)
#>  53°06.2869′ (S/W)
#>  25°14.4093′ (S/W)
#>   0°07.6754′ (S/W)
#>  18°23.9137′ (N/E)
#> 122°46.7203′ (S/W)
#>  77°02.1145′ (S/W)
#>   0°00.0000′ (N/E)
#>  12°17.3178′ (S/W)
#>  73°31.0370′ (N/E)
#>  57°31.1536′ (S/W)

## Name the "coords" object
names(dm) <- rep(c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
                   "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"), 2)
dm
#> Nelson's Column                           51°30.4659′ (N/E)
#> Ostravice                                 49°32.7726′ (N/E)
#> Tally Ho                                  48°06.4339′ (N/E)
#> Washington Monument                       38°53.3696′ (N/E)
#> Null Island                                0°00.0000′ (N/E)
#> Tristan da Cunha                          37°06.7044′ (S/W)
#> Mawson Peak                               53°06.2869′ (S/W)
#> Silvio Pettirossi International Airport   25°14.4093′ (S/W)
#> Nelson's Column                            0°07.6754′ (S/W)
#> Ostravice                                 18°23.9137′ (N/E)
#> Tally Ho                                 122°46.7203′ (S/W)
#> Washington Monument                       77°02.1145′ (S/W)
#> Null Island                                0°00.0000′ (N/E)
#> Tristan da Cunha                          12°17.3178′ (S/W)
#> Mawson Peak                               73°31.0370′ (N/E)
#> Silvio Pettirossi International Airport   57°31.1536′ (S/W)

## Set all values to represent longitude
## ("latlon" attribute set to FALSE, length 1)
latlon(dm) <- FALSE
dm
#> Nelson's Column                           51°30.4659′ E
#> Ostravice                                 49°32.7726′ E
#> Tally Ho                                  48°06.4339′ E
#> Washington Monument                       38°53.3696′ E
#> Null Island                                0°00.0000′ E
#> Tristan da Cunha                          37°06.7044′ W
#> Mawson Peak                               53°06.2869′ W
#> Silvio Pettirossi International Airport   25°14.4093′ W
#> Nelson's Column                            0°07.6754′ W
#> Ostravice                                 18°23.9137′ E
#> Tally Ho                                 122°46.7203′ W
#> Washington Monument                       77°02.1145′ W
#> Null Island                                0°00.0000′ E
#> Tristan da Cunha                          12°17.3178′ W
#> Mawson Peak                               73°31.0370′ E
#> Silvio Pettirossi International Airport   57°31.1536′ W

## Set eight values each of latitude and longitude
## ("latlon" attribute set to TRUE, n=8, and FALSE, n=8)
latlon(dm) <- rep(c(TRUE, FALSE), each = 8)
dm
#> Nelson's Column                           51°30.4659′ N
#> Ostravice                                 49°32.7726′ N
#> Tally Ho                                  48°06.4339′ N
#> Washington Monument                       38°53.3696′ N
#> Null Island                                0°00.0000′ N
#> Tristan da Cunha                          37°06.7044′ S
#> Mawson Peak                               53°06.2869′ S
#> Silvio Pettirossi International Airport   25°14.4093′ S
#> Nelson's Column                            0°07.6754′ W
#> Ostravice                                 18°23.9137′ E
#> Tally Ho                                 122°46.7203′ W
#> Washington Monument                       77°02.1145′ W
#> Null Island                                0°00.0000′ E
#> Tristan da Cunha                          12°17.3178′ W
#> Mawson Peak                               73°31.0370′ E
#> Silvio Pettirossi International Airport   57°31.1536′ W

## Show as an ordinary R numeric vector
as.numeric(dm)
#>  [1]   5130.4659   4932.7726   4806.4339   3853.3696      0.0000  -3706.7044
#>  [7]  -5306.2869  -2514.4093     -7.6754   1823.9137 -12246.7203  -7702.1145
#> [13]      0.0000  -1217.3178   7331.0370  -5731.1536

rm(dm)

```
