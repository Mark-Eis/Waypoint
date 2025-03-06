# Waypoint
### Conversion, Validation and Print Formatting of Geographic Coordinates.

The **Waypoint R package** enables conversion, validation and neatly formatted printing of
geographic coordinates and waypoints. Coordinates and waypoints are converted between (i) decimal
degrees, (ii) degrees and minutes, and (iii) degrees, minutes and seconds.

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

* Provides S3 `format()` and `print()` methods for neatly formatted printing of objects of classes
  "[`coords`](https://mark-eis.github.io/Waypoint/reference/coords.html)" and
  "[`waypoints`](https://mark-eis.github.io/Waypoint/reference/waypoints.html)".

*Waypoint* uses high performance C++ code seamlessly integrated into R using
[`Rcpp`](https://www.rcpp.org).
