# Extract or Replace Parts of a Coords Object

Extract or replace subsets of coords.

## Usage

``` r
# S3 method for class 'coords'
x[i]

# S3 method for class 'coords'
x[i] <- value
```

## Arguments

- x:

  a
  `"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`
  object.

- i:

  indices specifying elements to extract or replace—see base
  [`Extract`](https://rdrr.io/r/base/Extract.html).

- value:

  a `numeric`, a `numeric` vector of coordinate values of `length(i)`,
  or a `"coords"` object, possibly named.

## Value

A
`"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`
object.

## Details

Subsetting a
`"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`
object (except by an empty index) will drop all attributes except `fmt`,
`latlon`, `names` and `valid`. Indices referencing values greater than
`length(x)` will throw a `subscript out of bounds` error. If names are
not required, use [`unname()`](https://rdrr.io/r/base/unname.html), see
*examples*.

Replacement values may be a single `numeric`, a `numeric` vector of
coordinate values of `length(i)`, or a
`"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`
object, possibly with a `"latlon"` attribute. However, the `"latlon"`
attribute of the replacement value is ignored if the
`"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`
object `x` has no corresponding attribute set. If replacement values are
named, the names are also ignored; to replace names, use
[`names<-()`](https://rdrr.io/r/base/names.html) replacement form.

The `[<-(<coords>)` replacement operator automatically revalidates
`"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`
objects after the replacement operation by invoking
[`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md).

## Note

To extract and replace subsets of
`"`[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.md)`"`
objects, simply use the base package
[`[`](https://rdrr.io/r/base/Extract.html) and
[`[<-`](https://rdrr.io/r/base/Extract.html) operators, taking care not
to exclude the `latitude` and `longitude` columns or `"Name"` column (if
present), which could lead to undefined results. Whereas
`"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`
objects are automatically revalidated after using the `[<-(<coords>)`
replacement operator, following value replacement using the base
[`[<-`](https://rdrr.io/r/base/Extract.html) operator,
`"`[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.md)`"`
objects should be revalidated using
[`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md),
see
[`validate`](https://mark-eis.github.io/Waypoint/reference/validate.md)
*Examples*.

## See also

`"`[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.md)`"`,
base [`Extract`](https://rdrr.io/r/base/Extract.html),
[`unname()`](https://rdrr.io/r/base/unname.html), and
[`validate()`](https://mark-eis.github.io/Waypoint/reference/validate.md).

## Examples

``` r
## Continuing example from `as_coords()`...

## Named "coords" object in degrees and minutes with
## eight values each of latitude and longitude
dm
#> Nelson's Column                          51°30.4659′ N
#> Ostravice                                49°32.7726′ N
#> Tally Ho                                 48°06.4339′ N
#> Washington Monument                      38°53.3696′ N
#> Null Island                               0°00.0000′ N
#> Tristan da Cunha                         37°06.7044′ S
#> Mawson Peak                              53°06.2869′ S
#> Silvio Pettirossi International Airport  25°14.4093′ S
#> Nelson's Column                           0°07.6754′ W
#> Ostravice                                18°23.9137′ E
#> Tally Ho                                122°46.7203′ W
#> Washington Monument                      77°02.1145′ W
#> Null Island                               0°00.0000′ E
#> Tristan da Cunha                         12°17.3178′ W
#> Mawson Peak                              73°31.0370′ E
#> Silvio Pettirossi International Airport  57°31.1536′ W

## Extract the first eight values
dm[1:8]
#> Nelson's Column                          51°30.4659′ N
#> Ostravice                                49°32.7726′ N
#> Tally Ho                                 48°06.4339′ N
#> Washington Monument                      38°53.3696′ N
#> Null Island                               0°00.0000′ N
#> Tristan da Cunha                         37°06.7044′ S
#> Mawson Peak                              53°06.2869′ S
#> Silvio Pettirossi International Airport  25°14.4093′ S

## Exclude the first eight values
dm[-8:0]
#> Nelson's Column                           0°07.6754′ W
#> Ostravice                                18°23.9137′ E
#> Tally Ho                                122°46.7203′ W
#> Washington Monument                      77°02.1145′ W
#> Null Island                               0°00.0000′ E
#> Tristan da Cunha                         12°17.3178′ W
#> Mawson Peak                              73°31.0370′ E
#> Silvio Pettirossi International Airport  57°31.1536′ W

## Index odd-numbered values
(index <- as.logical(1:16 %% 2))
#>  [1]  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE
#> [13]  TRUE FALSE  TRUE FALSE
dm[index]
#> Nelson's Column  51°30.4659′ N
#> Tally Ho         48°06.4339′ N
#> Null Island       0°00.0000′ N
#> Mawson Peak      53°06.2869′ S
#> Nelson's Column   0°07.6754′ W
#> Tally Ho        122°46.7203′ W
#> Null Island       0°00.0000′ E
#> Mawson Peak      73°31.0370′ E

## Extract values without names
unname(dm)[1:4]
#>  51°30.4659′ N
#>  49°32.7726′ N
#>  48°06.4339′ N
#>  38°53.3696′ N

## Create "coords" object with updated position of Tally Ho
newpos <- as_coords(c(4930.342, -12411.580), fmt = 2)
latlon(newpos) <- c(TRUE, FALSE)
newpos
#>  49°30.3420′ N
#> 124°11.5800′ W

## Update position using the "coords" object as replacement value
dm[c(3, 11)] <- newpos
dm[c(3, 11)]
#> Tally Ho  49°30.3420′ N
#> Tally Ho 124°11.5800′ W

## Or, as latlon didn't actually change, use simple numeric vector
dm[c(3, 11)] <- c(4930.342, -12411.580)
dm[c(3, 11)]
#> Tally Ho  49°30.3420′ N
#> Tally Ho 124°11.5800′ W

rm(dm, index, newpos)
```
