# Convert the Format of "coords" and "waypoints" Objects

Convert the format of objects of class `"coords"` or `"waypoints"`
between (i) decimal degrees, (ii) degrees and minutes, and (iii)
degrees, minutes and seconds.

## Usage

``` r
convert(x, ...)

# S3 method for class 'coords'
convert(x, fmt, ...)

# S3 method for class 'waypoints'
convert(x, fmt, ...)
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

- fmt:

  `integer`, `1L`, `2L` or `3L`, specifying the required coordinate
  format.

## Value

The original argument `x`, an object of class `"coords"` or
`"waypoints"` with values converted as described under *details* and a
revised `"fmt"` attribute reflecting the new format.

## Details

The `fmt` argument should be `1L` to convert to *decimal degrees*, `2L`,
to convert to *degrees and minutes*, and `3L` to convert to *degrees,
minutes and seconds*. On conversion of a `"coords"` object, the original
argument `x` is modified to have a decimal point after the number of
whole degrees in the case of decimal degrees, after the number of whole
minutes in the case of degrees and minutes, and after the number of
whole seconds in the case of degrees, minutes and seconds.

Prior to conversion, the `"coords"` or `"waypoints"` object to be
converted is checked to ensure its values represent valid geographic
locations as described under
[`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md).

## Note

`convert()` modifies its argument `x` in place. To format or print
`"coords"` or `"waypoints"` in another coordinate format without
modifying the original object, use
[`format()`](https://mark-eis.github.io/Waypoint/reference/format.md) or
[`print()`](https://mark-eis.github.io/Waypoint/reference/format.md).

## See also

`"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`,
`"`[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.md)`"`
and
[`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md).

Other coordsandway:
[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md),
[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.md)

## Examples

``` r
## Continuing example from `as_coords()`...

## Named "coords" object in degrees and minutes with
## eight values each of latitude and longitude
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

## Convert to degrees, minutes and seconds (fmt = 3)
convert(dm, 3)
#> Nelson's Column                           51°30′27.95″ N
#> Ostravice                                 49°32′46.36″ N
#> Tally Ho                                  48°06′26.03″ N
#> Washington Monument                       38°53′22.18″ N
#> Null Island                                0°00′00.00″ N
#> Tristan da Cunha                          37°06′42.26″ S
#> Mawson Peak                               53°06′17.21″ S
#> Silvio Pettirossi International Airport   25°14′24.56″ S
#> Nelson's Column                            0°07′40.52″ W
#> Ostravice                                 18°23′54.82″ E
#> Tally Ho                                 122°46′43.22″ W
#> Washington Monument                       77°02′06.87″ W
#> Null Island                                0°00′00.00″ E
#> Tristan da Cunha                          12°17′19.07″ W
#> Mawson Peak                               73°31′02.22″ E
#> Silvio Pettirossi International Airport   57°31′09.22″ W

## Convert to decimal degrees (fmt = 1)
convert(dm, 1)
#> Nelson's Column                            51.507765° lat
#> Ostravice                                  49.546210° lat
#> Tally Ho                                   48.107232° lat
#> Washington Monument                        38.889493° lat
#> Null Island                                 0.000000° lat
#> Tristan da Cunha                          -37.111740° lat
#> Mawson Peak                               -53.104782° lat
#> Silvio Pettirossi International Airport   -25.240155° lat
#> Nelson's Column                            -0.127923° lon
#> Ostravice                                  18.398562° lon
#> Tally Ho                                 -122.778672° lon
#> Washington Monument                       -77.035242° lon
#> Null Island                                 0.000000° lon
#> Tristan da Cunha                          -12.288630° lon
#> Mawson Peak                                73.517283° lon
#> Silvio Pettirossi International Airport   -57.519227° lon

## Show converted values as an ordinary R numeric vector
as.numeric(dm)
#>  [1]   51.5077650   49.5462100   48.1072317   38.8894933    0.0000000
#>  [6]  -37.1117400  -53.1047817  -25.2401550   -0.1279233   18.3985617
#> [11] -122.7786717  -77.0352417    0.0000000  -12.2886300   73.5172833
#> [16]  -57.5192267

###
## Continuing example from `as_waypoints()`...

## "waypoints" object in degrees, minutes and seconds
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

## Convert to degrees and minutes (fmt = 2)
convert(wp, 2)
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

## Convert to decimal degrees (fmt = 1)
convert(wp, 1)
#>                                             Latitude     Longitude
#>                                           ___________  ____________
#> Nelson's Column                            51.507764°    -0.127925°
#> Ostravice                                  49.546211°    18.398561°
#> Tally Ho                                   48.107233°  -122.778672°
#> Washington Monument                        38.889494°   -77.035242°
#> Null Island                                 0.000000°     0.000000°
#> Tristan da Cunha                          -37.111739°   -12.288631°
#> Mawson Peak                               -53.104781°    73.517283°
#> Silvio Pettirossi International Airport   -25.240156°   -57.519225°

## Show converted values as an ordinary R data frame
as.data.frame(wp)
#>                                      name       lat         lon
#> 1                         Nelson's Column  51.50776   -0.127925
#> 2                               Ostravice  49.54621   18.398561
#> 3                                Tally Ho  48.10723 -122.778672
#> 4                     Washington Monument  38.88949  -77.035242
#> 5                             Null Island   0.00000    0.000000
#> 6                        Tristan da Cunha -37.11174  -12.288631
#> 7                             Mawson Peak -53.10478   73.517283
#> 8 Silvio Pettirossi International Airport -25.24016  -57.519225

rm(dm, wp)

```
