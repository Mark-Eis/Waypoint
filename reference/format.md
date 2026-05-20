# Format and Print Coords or Waypoints

Format and print objects of class `"coords"` or `"waypoints"`.

## Usage

``` r
# S3 method for class 'coords'
print(x, ..., max = NULL)

# S3 method for class 'waypoints'
print(x, ..., fmt = NULL, max = NULL)

# S3 method for class 'coords'
format(x, ..., usenames = TRUE, validate = TRUE, fmt = 0L)

# S3 method for class 'waypoints'
format(x, ..., usenames = TRUE, validate = TRUE, fmt = 0L)

ll_headers(width, fmt)
```

## Arguments

- x:

  object of class
  `"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`
  created by function
  [`as_coords()`](https://mark-eis.github.io/Waypoint/reference/coords.md),
  or class
  `"`[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.md)`"`
  created by function
  [`as_waypoints()`](https://mark-eis.github.io/Waypoint/reference/waypoints.md).

- ...:

  further arguments passed to or from other methods.

- max:

  numeric or `NULL`, specifying the maximal number of entries to be
  printed. By default, when `NULL`,
  [`getOption`](https://rdrr.io/r/base/options.html)`("max.print")`
  used.

- fmt:

  `integer`, `1L`, `2L` or `3L`, specifying the required coordinate
  format.

- usenames:

  `logical`, whether or not to include names in formatted output;
  default `TRUE`.

- validate:

  `logical`, whether or not to
  [`validate`](https://mark-eis.github.io/Waypoint/reference/validate.md)
  `x` before formatting; default `TRUE`.

- width:

  `character` vector, used to match width of headers to formatted
  output.

## Value

The `format()` methods for both classes `"coords"` and `"waypoints"`
return a `character` vector, respectively of length `length(x)` or
`nrow(x)`, and containing values formatted in decimal degrees, degrees
and minutes, or degrees, minutes and seconds as appropriate.

## Details

The `format()` methods for
`"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`
and
`"`[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.md)`"`
objects output elegantly formatted `character` vector representations of
their arguments, which are used by their respective
[`print()`](https://rdrr.io/r/base/print.html) methods.

Objects of class `"coords"` specified in *degrees and minutes* or in
*degrees, minutes and seconds* and with a `"latlon"` attribute, and
similarly specified `"waypoints"` objects are formatted with individual
coordinate values followed by a capital letter representing the
*cardinal direction* i.e., `N`, `E`, `S` or `W`. `"coords"` objects
lacking a `"latlon"` attribute have formatted values followed by two
possible cardinal directions in parentheses i.e., `(N/E)` for positive
values and `(S/W)` for negative values. Values of `"coords"` or
`"waypoints"` objects in *decimal degrees* are formatted prefixed with
their sign, if negative; cardinal direction is not shown, but for
`"coords"` objects with a `"latlon"` attribute, the formatted values are
suffixed by either `lat` or `lon`.

Prior to formatting and printing, `"coords"` or `"waypoints"` objects
are checked to ensure that their `"valid"` attribute (in the case of a
`"coords"` object), or `"validlat"` and `"validlon"` attributes (in the
case of a `"waypoints"` object) are present and all `TRUE` i.e., valid.
If these attributes are found to contain any `FALSE` i.e. invalid
values, a warning is issued and similarly, if these attributes are
missing, a warning is issued and the objects are re-validated as
described under
[`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md).

The optional argument `fmt` may be used to specify the coordinate format
desired for formatting or printing `"coords"` or `"waypoints"` objects,
see the `fmt` argument for
[`as_coords()`](https://mark-eis.github.io/Waypoint/reference/coords.md)
and
[`as_waypoints()`](https://mark-eis.github.io/Waypoint/reference/waypoints.md);
using the default, `fmt = 0L`, will format or print in the existing
coordinate format.

`ll_headers()` outputs the headings `"Latitude ... Longitude"` formatted
to the width of argument `width`, adjusted for format `fmt` and is
primarily intended for use by the
[`print()`](https://rdrr.io/r/base/print.html) method for class
`"waypoints"`. Likewise argument `validate` is used by the
[`print()`](https://rdrr.io/r/base/print.html) methods for classes
`"coords"` and `"waypoints"` to prevent unnecessary replicate validation
and may otherwise be left as the default.

## See also

[`format()`](https://rdrr.io/r/base/format.html),
[`print()`](https://rdrr.io/r/base/print.html),
`"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`
and
`"`[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.md)`"`.

## Examples

``` r
## Continuing example from `as_coords()`...

## Print named "coords" object in degrees and minutes,
## implicitly using S3 print() method
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

## Print explicitly using S3 print() method, specifying
## the maximal number of entries to be printed
print(dm, max = 14)
#> Nelson's Column       51°30.4659′ N
#> Ostravice             49°32.7726′ N
#> Tally Ho              48°06.4339′ N
#> Washington Monument   38°53.3696′ N
#> Null Island            0°00.0000′ N
#> Tristan da Cunha      37°06.7044′ S
#> Mawson Peak           53°06.2869′ S
#>  [ reached 'max' / getOption("max.print") -- omitted 9 entries ]

## Format as a fixed-width character vector,
## with names...
format(dm)
#>  [1] "Nelson's Column                           51°30.4659′ N"
#>  [2] "Ostravice                                 49°32.7726′ N"
#>  [3] "Tally Ho                                  48°06.4339′ N"
#>  [4] "Washington Monument                       38°53.3696′ N"
#>  [5] "Null Island                                0°00.0000′ N"
#>  [6] "Tristan da Cunha                          37°06.7044′ S"
#>  [7] "Mawson Peak                               53°06.2869′ S"
#>  [8] "Silvio Pettirossi International Airport   25°14.4093′ S"
#>  [9] "Nelson's Column                            0°07.6754′ W"
#> [10] "Ostravice                                 18°23.9137′ E"
#> [11] "Tally Ho                                 122°46.7203′ W"
#> [12] "Washington Monument                       77°02.1145′ W"
#> [13] "Null Island                                0°00.0000′ E"
#> [14] "Tristan da Cunha                          12°17.3178′ W"
#> [15] "Mawson Peak                               73°31.0370′ E"
#> [16] "Silvio Pettirossi International Airport   57°31.1536′ W"

## ...or without them
format(dm, usenames = FALSE)
#>  [1] " 51°30.4659′ N" " 49°32.7726′ N" " 48°06.4339′ N" " 38°53.3696′ N"
#>  [5] "  0°00.0000′ N" " 37°06.7044′ S" " 53°06.2869′ S" " 25°14.4093′ S"
#>  [9] "  0°07.6754′ W" " 18°23.9137′ E" "122°46.7203′ W" " 77°02.1145′ W"
#> [13] "  0°00.0000′ E" " 12°17.3178′ W" " 73°31.0370′ E" " 57°31.1536′ W"

## Format as decimal degrees,
format(dm, fmt = 1)
#>  [1] "Nelson's Column                            51.507765° lat"
#>  [2] "Ostravice                                  49.546210° lat"
#>  [3] "Tally Ho                                   48.107232° lat"
#>  [4] "Washington Monument                        38.889493° lat"
#>  [5] "Null Island                                 0.000000° lat"
#>  [6] "Tristan da Cunha                          -37.111740° lat"
#>  [7] "Mawson Peak                               -53.104782° lat"
#>  [8] "Silvio Pettirossi International Airport   -25.240155° lat"
#>  [9] "Nelson's Column                            -0.127923° lon"
#> [10] "Ostravice                                  18.398562° lon"
#> [11] "Tally Ho                                 -122.778672° lon"
#> [12] "Washington Monument                       -77.035242° lon"
#> [13] "Null Island                                 0.000000° lon"
#> [14] "Tristan da Cunha                          -12.288630° lon"
#> [15] "Mawson Peak                                73.517283° lon"
#> [16] "Silvio Pettirossi International Airport   -57.519227° lon"

###
## Continuing example from `as_waypoints()`...

## Print named "waypoints" object in degrees, minutes and seconds
## implicitly using S3 print() method
wp
#>                                                Latitude        Longitude
#>                                           ______________  _______________
#> Nelson's Column                           51°30′27.95″ N    0°07′40.53″ W
#> Ostravice                                 49°32′46.36″ N   18°23′54.82″ E
#> Tally Ho                                  48°06′26.04″ N  122°46′43.22″ W
#> Washington Monument                       38°53′22.18″ N   77°02′06.87″ W
#> Null Island                                0°00′00.00″ N    0°00′00.00″ E
#> Tristan da Cunha                          37°06′42.26″ S   12°17′19.07″ W
#> Mawson Peak                               53°06′17.21″ S   73°31′02.22″ E
#> Silvio Pettirossi International Airport   25°14′24.56″ S   57°31′09.21″ W

## Print explicitly using S3 print() method, specifying
## the maximal number of entries to be printed
print(wp, max = 21)
#>                            Latitude        Longitude
#>                       ______________  _______________
#> Nelson's Column       51°30′27.95″ N    0°07′40.53″ W
#> Ostravice             49°32′46.36″ N   18°23′54.82″ E
#> Tally Ho              48°06′26.04″ N  122°46′43.22″ W
#> Washington Monument   38°53′22.18″ N   77°02′06.87″ W
#> Null Island            0°00′00.00″ N    0°00′00.00″ E
#> Tristan da Cunha      37°06′42.26″ S   12°17′19.07″ W
#> Mawson Peak           53°06′17.21″ S   73°31′02.22″ E
#>  [ reached 'max' / getOption("max.print") -- omitted 1 rows ]

## Print as degrees and minutes
print(wp, fmt = 2)
#>                                               Latitude       Longitude
#>                                           _____________  ______________
#> Nelson's Column                           51°30.4658′ N    0°07.6755′ W
#> Ostravice                                 49°32.7727′ N   18°23.9137′ E
#> Tally Ho                                  48°06.4340′ N  122°46.7203′ W
#> Washington Monument                       38°53.3697′ N   77°02.1145′ W
#> Null Island                                0°00.0000′ N    0°00.0000′ E
#> Tristan da Cunha                          37°06.7043′ S   12°17.3178′ W
#> Mawson Peak                               53°06.2868′ S   73°31.0370′ E
#> Silvio Pettirossi International Airport   25°14.4093′ S   57°31.1535′ W

## Format as a fixed-width character vector,
## with names...
format(wp)
#> [1] "Nelson's Column                           51°30′27.95″ N    0°07′40.53″ W"
#> [2] "Ostravice                                 49°32′46.36″ N   18°23′54.82″ E"
#> [3] "Tally Ho                                  48°06′26.04″ N  122°46′43.22″ W"
#> [4] "Washington Monument                       38°53′22.18″ N   77°02′06.87″ W"
#> [5] "Null Island                                0°00′00.00″ N    0°00′00.00″ E"
#> [6] "Tristan da Cunha                          37°06′42.26″ S   12°17′19.07″ W"
#> [7] "Mawson Peak                               53°06′17.21″ S   73°31′02.22″ E"
#> [8] "Silvio Pettirossi International Airport   25°14′24.56″ S   57°31′09.21″ W"

## ...or without them
format(wp, usenames = FALSE)
#> [1] " 51°30′27.95″ N    0°07′40.53″ W" " 49°32′46.36″ N   18°23′54.82″ E"
#> [3] " 48°06′26.04″ N  122°46′43.22″ W" " 38°53′22.18″ N   77°02′06.87″ W"
#> [5] "  0°00′00.00″ N    0°00′00.00″ E" " 37°06′42.26″ S   12°17′19.07″ W"
#> [7] " 53°06′17.21″ S   73°31′02.22″ E" " 25°14′24.56″ S   57°31′09.21″ W"

rm(dm, wp)

```
