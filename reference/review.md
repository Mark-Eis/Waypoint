# Review Coordinates and Waypoints Validity

Review validity of elements of `"coords"` and `"waypoints"` objects.

## Usage

``` r
review(x, ...)

# S3 method for class 'coords'
review(x, ..., show_n = 20L)

# S3 method for class 'waypoints'
review(x, ..., show_n = 20L)
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

- show_n:

  `integer`, the maximum number of invalid elements of argument `x` to
  include in the output; default `20L`.

## Value

The `review()` method for class `"coords"` returns a
[`list`](https://rdrr.io/r/base/list.html) comprising the following
elements: -

- allvalid:

  `logical`, whether or not all the elements of argument `x` are valid.

- n_invalid:

  `integer`, the number of invalid elements in argument `x`, if any.

- invalids:

  `numeric` vector including invalid elements of argument `x`, if any.

- which_invalid:

  `integer` vector specifying which elements of argument `x` are
  invalid, if any.

The method for class `"waypoints"` returns a list of two sub-lists, each
sub-list with elements as described above for the method for class
`"coords"`, one each for latitude and longitude.

## Details

`review()` reveals elements of `"coords"` and `"waypoints"` objects that
do not conform to the criteria checked by
[`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md),
i.e. are not valid geographic locations.

## See also

`"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`
and
`"`[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.md)`"`.

Other validate:
[`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md)

## Examples

``` r
## Continuing example from `validate()`...

## Review "coords" object in degrees and minutes, having
## an erroneous first value of more than 60 minutes
review(dm)
#> $allvalid
#> [1] FALSE
#> 
#> $n_invalid
#> [1] 1
#> 
#> $invalids
#> Warning: Invalid coords!
#> Nelson's Column   51°60.4659′ N
#> 
#> $which_invalid
#> [1] 1
#> 

###
## Continuing example from `validate()`...

## Review "waypoints" object in  decimal degrees, having an erroneous
## penultimate latitude absolute value greater than 90 degrees
review(wp)
#> $Lat
#> $Lat$allvalid
#> [1] FALSE
#> 
#> $Lat$n_invalid
#> [1] 1
#> 
#> $Lat$invalids
#> Warning: Invalid coords!
#> Mawson Peak   -93.104781° lat
#> 
#> $Lat$which_invalid
#> [1] 7
#> 
#> 
#> $Lon
#> $Lon$allvalid
#> [1] TRUE
#> 
#> $Lon$n_invalid
#> [1] 0
#> 
#> $Lon$invalids
#> [1] NA
#> 
#> $Lon$which_invalid
#> [1] NA
#> 
#> 

rm(dm, wp)
```
