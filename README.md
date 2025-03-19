# Waypoint
### Conversion, Validation, Formatting and Printing of Geographic Coordinates

The **Waypoint R package** enables conversion, validation, formatting and and elegant printing of
geographic positional coordinates and waypoints (paired latitude and longitude values).
Coordinates and waypoints may be converted easily and rapidly between (i) decimal degree, (ii)
degrees and minutes, and (iii) degrees, minutes and seconds formats.

## Installation

You can install the currently-released version from CRAN with this R
command:

``` r
install.packages("Waypoint")
```

Alternatively, you can install the latest development version of Waypoint
from [GitHub](https://github.com/) with:
      
``` r
# install.packages("devtools")
devtools::install_github("Mark-Eis/Waypoint")
```
---

**Author:** Mark C. Eisler

**eMail:** Mark.Eisler@bristol.ac.uk

**ORCID** = [0000-0001-6843-3345](https://orcid.org/0000-0001-6843-3345)

### Waypoint Package Overview: –

* Create "[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.html)"
  objects in each format with `as_coords()`.

* Create
 "[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.html)"
 objects in each format with `as_waypoints()`.

* Convert "`coords`" and "`waypoints`" objects between decimal degrees, degrees
  and minutes, and degrees, minutes and seconds formats with `convert()`.

* Assign latitude and longitude attributes to individual coordinate values
  within "`coords`" objects with
 [`latlon<-`](https://mark-eis.github.io/Waypoint/reference/latlon.html)().

* Ensure values within "`coords`" and "`waypoints`" objects are valid
  geographic locations with `validate()` and identifies individual invalid
  values with `review()`.

* Use `format()` and `print()` S3 methods for neat formatting and printing
  of objects of classes "`coords`" and "`waypoints`".
  
#### Methodology  

*Waypoint* uses high performance C++ code seamlessly integrated into R using
[`Rcpp`](https://www.rcpp.org) to enable rapid conversion and formatting of
large coordinate and waypoint datasets.

#### Disclaimer

While every effort is made to ensure this package functions as expected, the
author accepts no responsibility for the consequences of errors if your map
shows a city in the middle of the ocean, if your boat runs aground, or even
if your aeroplane crashes into the mountain.
