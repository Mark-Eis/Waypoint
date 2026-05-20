# Validate Coords or Waypoints

Validate objects of class `"coords"` or `"waypoints"` as geographic
locations.

## Usage

``` r
validate(x, ...)

# S3 method for class 'coords'
validate(x, ..., force = TRUE)

# S3 method for class 'waypoints'
validate(x, ..., force = TRUE)
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

- force:

  `logical` signifying whether, if `TRUE`, to perform full *de novo*
  revalidation or, if `FALSE`, simply check existing `"valid"` attribute
  in the case of a `"coords"` object, or `"validlat"` and `"validlon"`
  attributes in the case of a `"waypoints"` object and only revalidate
  if any of these are missing; default `TRUE`.

## Value

`validate()` returns its argument with `logical` vector attribute
`"valid"`, or attributes `"validlat"` and `"validlon"` updated as
appropriate for `"coords"` and' `"waypoints"` objects respectively.

## Details

Individual coordinate values within
`"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`
or
`"`[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.md)`"`
objects are checked to ensure they represent valid geographic locations.

To be valid, the absolute values of coordinates in degrees must not
exceed 180°, or 90° if degrees of latitude and, similarly, the absolute
values of the minutes and seconds components, where given, must not
exceed 60. Otherwise, a warning will be issued and the `"valid"`
attribute in the case of a `"coords"` object, or `"validlat"` and
`"validlon"` attributes in the case of a `"waypoints"` object will be
set to `FALSE` for any non-compliant coordinate values.

Argument `force` is primarily intended for use by the
[`print()`](https://rdrr.io/r/base/print.html) methods for classes
`"coords"` and `"waypoints"` and should otherwise left as the default
value `TRUE`.

## See also

`"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`
and
`"`[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.md)`"`.

Other validate:
[`review()`](https://mark-eis.github.io/Waypoint/reference/review.md)

## Examples

``` r
## Continuing example from `as_coords()`...

## Validate "coords" object in degrees and minutes
validate(dm)
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

## Deliberately change the first coordinate
## to a value greater than 60 minutes
dm[1] <- 5160.4659
#> Warning: Validation failed!

validate(dm)
#> Warning: Validation failed!
#> Warning: Invalid coords!
#> Nelson's Column                           51°60.4659′ N
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

## Examine "valid" attribute of dm
attr(dm, "valid")
#>  [1] FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#> [13]  TRUE  TRUE  TRUE  TRUE

###
## Continuing second example from `as_waypoints()`...

## Validate "waypoints" object in decimal degrees

validate(wp)
#>                                             Latitude     Longitude
#>                                           ___________  ____________
#> Nelson's Column                            51.507765°    -0.127924°
#> Ostravice                                  49.546210°    18.398562°
#> Tally Ho                                   48.107232°  -122.778671°
#> Washington Monument                        38.889494°   -77.035242°
#> Null Island                                 0.000000°     0.000000°
#> Tristan da Cunha                          -37.111740°   -12.288630°
#> Mawson Peak                               -53.104781°    73.517283°
#> Silvio Pettirossi International Airport   -25.240156°   -57.519227°

## Deliberately change the penultimate latitude
## to an absolute value greater than 90 degrees
wp$lat[7] <- -93.104781

validate(wp)
#> Warning: Validation of latitude failed!
#> Warning: Invalid latitude!
#>                                             Latitude     Longitude
#>                                           ___________  ____________
#> Nelson's Column                            51.507765°    -0.127924°
#> Ostravice                                  49.546210°    18.398562°
#> Tally Ho                                   48.107232°  -122.778671°
#> Washington Monument                        38.889494°   -77.035242°
#> Null Island                                 0.000000°     0.000000°
#> Tristan da Cunha                          -37.111740°   -12.288630°
#> Mawson Peak                               -93.104781°    73.517283°
#> Silvio Pettirossi International Airport   -25.240156°   -57.519227°

## Examine "validlat" attribute of wp
attr(wp, "validlat")
#> [1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE

rm(dm, wp)
```
