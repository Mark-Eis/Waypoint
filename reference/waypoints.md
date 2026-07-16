# Geographic Waypoint Class

`as_waypoints()` creates an object of class `"waypoints"`, a robust
representation of a series of geographic or GPS waypoints of paired
latitude and longitude values.

## Usage

``` r
as_waypoints(object, ...)

# Default S3 method
as_waypoints(object, ..., fmt = 1L)
```

## Arguments

- object:

  a data frame with each row representing a waypoint, comprising at
  least two `numeric` columns containing values of latitude and
  longitude, and optionally a `character` column of waypoint names (see
  *Details*).

- ...:

  further arguments passed to or from other methods.

- fmt:

  `integer`, `1L`, `2L` or `3L`, specifying the required coordinate
  format.

## Value

An object of classes `"waypoints"` and `"data.frame"`, comprising the
original data frame argument `object`, with additional attributes: –

- `"class"`:

  the `character` string "waypoints".

- `"fmt"`:

  an `integer` indicating the coordinate format.

- `"namescol"`:

  an `integer` indicating the position of a waypoint names column, if
  present.

- `"llcols"`:

  a length 2 `integer` vector indicating the positions of latitude and
  longitude columns.

- `"validlat"` and `"validlon"`:

  `logical` vectors indicating whether individual latitude and longitude
  values are valid geographic locations.

## Details

Data frame argument `object` should have `numeric` vector latitude and
longitude columns with individual values having a decimal point after
the number of whole degrees in the case of *decimal degrees*, after the
number of whole minutes in the case of *degrees and minutes*, and after
the number of whole seconds in the case of *degrees, minutes and
seconds*. These should be the first two columns of the data frame, or
the second and third columns if the first column contains waypoints
names (see below). Alternative columns may be specified for the latitude
and longitude by setting `"llcols"` as a length 2 `integer` vector
attribute of `object` indicating their positions in the data frame.

The `fmt` argument should be `1L` to represent decimal degrees, `2L` for
degrees and minutes, and `3L` for degrees, minutes and seconds and is
used to provide the format of values in data frame argument `object` to
be converted to class `"waypoints"`.

If the waypoints have names, these should be included in a "Name" column
of data frame argument `object`, by default immediately before (on the
left-hand side of) the latitude and longitude columns. An alternative
column for waypoint names may be specified by setting an `integer`
[`attribute`](https://rdrr.io/r/base/attributes.html) named "namescol"
indicating its position in `object`. If neither a "Name" column nor a
`"namescol"` attribute is present in `object`,
[`row.names`](https://rdrr.io/r/base/row.names.html) are used for
waypoint names.

The latitude and longitude values of a newly created `"waypoints"`
object are checked to ensure they are valid geographic locations as
described under
[`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md).

## See also

[`attr()`](https://rdrr.io/r/base/attr.html),
[`data.frame()`](https://rdrr.io/r/base/data.frame.html)
[`format`](https://mark-eis.github.io/Waypoint/reference/format.md)`(<waypoints>)`,
[`print`](https://mark-eis.github.io/Waypoint/reference/format.md)`(<waypoints>)`,
and
[`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md).

Other coordsandway:
[`convert()`](https://mark-eis.github.io/Waypoint/reference/convert.md),
[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)

## Examples

``` r
## Dataframe representing waypoint names, and latitude and longitude values
## in degrees, minutes and seconds
wp1 <- data.frame(
    name = c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
             "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"),
    lat = c(513027.95, 493246.36, 480626.04, 385322.18, 0, -370642.26, -530617.21, -251424.56),
    lon = c(-00740.53, 182354.82, -1224643.22, -770206.87, 0, -121719.07, 733102.22, -573109.21)
)

## Create "waypoints" object in degrees, minutes and seconds (fmt = 3)
as_waypoints(wp1, fmt = 3)
#>                                               Latitude        Longitude
#>                                          ______________  _______________
#> Nelson's Column                          51°30′27.95″ N    0°07′40.53″ W
#> Ostravice                                49°32′46.36″ N   18°23′54.82″ E
#> Tally Ho                                 48°06′26.04″ N  122°46′43.22″ W
#> Washington Monument                      38°53′22.18″ N   77°02′06.87″ W
#> Null Island                               0°00′00.00″ N    0°00′00.00″ E
#> Tristan da Cunha                         37°06′42.26″ S   12°17′19.07″ W
#> Mawson Peak                              53°06′17.21″ S   73°31′02.22″ E
#> Silvio Pettirossi International Airport  25°14′24.56″ S   57°31′09.21″ W

## Show as an ordinary R data frame
as.data.frame(wp1)
#>                                      name       lat         lon
#> 1                         Nelson's Column  513028.0     -740.53
#> 2                               Ostravice  493246.4   182354.82
#> 3                                Tally Ho  480626.0 -1224643.22
#> 4                     Washington Monument  385322.2  -770206.87
#> 5                             Null Island       0.0        0.00
#> 6                        Tristan da Cunha -370642.3  -121719.07
#> 7                             Mawson Peak -530617.2   733102.22
#> 8 Silvio Pettirossi International Airport -251424.6  -573109.21

###

## Dataframe representing unnamed latitude and longitude
## values in decimal degrees
wp2 <- data.frame(
    lat = c(51.507765, 49.54621, 48.107232, 38.889494, 0, -37.11174, -53.104781, -25.240156),
    lon = c(-0.127924, 18.398562, -122.778671, -77.035242, 0, -12.28863, 73.517283, -57.519227)
)

## Create unnamed "waypoints" object in decimal degrees (default fmt = 1)
as_waypoints(wp2)
#>      Latitude     Longitude
#>    ___________  ____________
#> 1   51.507765°    -0.127924°
#> 2   49.546210°    18.398562°
#> 3   48.107232°  -122.778671°
#> 4   38.889494°   -77.035242°
#> 5    0.000000°     0.000000°
#> 6  -37.111740°   -12.288630°
#> 7  -53.104781°    73.517283°
#> 8  -25.240156°   -57.519227°

## Add waypoint names as row.names
row.names(wp2) <-
    c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
      "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport")

wp2
#>                                            Latitude     Longitude
#>                                          ___________  ____________
#> Nelson's Column                           51.507765°    -0.127924°
#> Ostravice                                 49.546210°    18.398562°
#> Tally Ho                                  48.107232°  -122.778671°
#> Washington Monument                       38.889494°   -77.035242°
#> Null Island                                0.000000°     0.000000°
#> Tristan da Cunha                         -37.111740°   -12.288630°
#> Mawson Peak                              -53.104781°    73.517283°
#> Silvio Pettirossi International Airport  -25.240156°   -57.519227°

## Show as an ordinary R data frame
as.data.frame(wp2)
#>                                               lat         lon
#> Nelson's Column                          51.50776   -0.127924
#> Ostravice                                49.54621   18.398562
#> Tally Ho                                 48.10723 -122.778671
#> Washington Monument                      38.88949  -77.035242
#> Null Island                               0.00000    0.000000
#> Tristan da Cunha                        -37.11174  -12.288630
#> Mawson Peak                             -53.10478   73.517283
#> Silvio Pettirossi International Airport -25.24016  -57.519227

rm(wp1, wp2)

```
