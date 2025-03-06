# Waypoint
### Conversion, Validation, Formatting and Printing of Geographic Coordinates

The **Waypoint R package** enables conversion, validation, and neat formatting and printing of
geographic positional coordinate values and waypoints of paired latitude and longitude values.
Coordinates and waypoints may be converted easily and rapidly between (i) decimal degree, (ii)
degrees and minutes, and (iii) degrees, minutes and seconds formats.

## Installation

You can install the development version of Waypoint from [GitHub](https://github.com/) with:
      
``` r
# install.packages("devtools")
devtools::install_github("Mark-Eis/Waypoint")
```
---

**Author:** Mark C. Eisler

**eMail:** Mark.Eisler@bristol.ac.uk

**ORCID** = [0000-0001-6843-3345](https://orcid.org/0000-0001-6843-3345)

### Waypoint Package Overview: –

* Creates "[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.html)" objects in each
  format with `as_coords`().

* Creates "[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.html)" objects in
  each format with `as_waypoints`().

* Converts "[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.html)" and
  "[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.html)" objects
  between decimal degrees, degrees and minutes, and degrees, minutes and seconds formats with
  `convert()`.

* Assigns latitude and longitude attributes to individual coordinate values within
  "[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.html)" objects with
  [`latlon<-`](https://mark-eis.github.io/Waypoint/reference/latlon.html)().

* Ensures values within "[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.html)" and
  "[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.html)" objects are valid
  geographic locations with `validate()` and identifies individual invalid values with `review()`.

* Provides `format()` and `print()` S3 methods for neat formatting and printing of objects of
  classes "[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.html)" and
  "[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.html)".
  
#### Methodology  

*Waypoint* uses high performance C++ code seamlessly integrated into R using
[`Rcpp`](https://www.rcpp.org) to enable rapid conversion and formatting of large coordinate and
waypoint datasets.

#### Disclaimer

While every effort is made to ensure this package functions as expected, the author accepts no
responsibility for the consequences of errors even if your map shows your city in the middle of the
ocean, your boat runs aground, or your aeroplane crashes into the mountain.
