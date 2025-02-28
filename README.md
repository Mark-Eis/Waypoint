# Waypoint
### Conversion, Validation and Print Formatting of Geographic Coordinates.

**Author:** Mark C. Eisler

**eMail:** Mark.Eisler@bristol.ac.uk

**ORCID** = [0000-0001-6843-3345](https://orcid.org/0000-0001-6843-3345)

## Installation

You can install the development version of Waypoint from [GitHub](https://github.com/) with:
      
``` r
# install.packages("devtools")
devtools::install_github("Mark-Eis/Waypoint")
```
---

### Waypoint Package Description: –

The **Waypoint R package** enables conversion, validation and neatly formatted printing of geographic coordinates and waypoints. Coordinates and waypoints are converted between (i) decimal degrees, (ii) degrees and minutes, and (iii) degrees, minutes and seconds. More specifically, *Waypoint* does the following: – 

- Creates [`"coords"`](https://mark-eis.github.io/Waypoint/reference/coords.html) objects and converts their formats with `coords()` and [`coords<-()`](https://mark-eis.github.io/Waypoint/reference/coords.html).

- Creates [`"waypoints"`](https://mark-eis.github.io/Waypoint/reference/waypoints.html) objects and converts their formats with `waypoints()` and  [`waypoints<-()`](https://mark-eis.github.io/Waypoint/reference/waypoints.html).

- Assigns latitude and longitude attributes to values within objects of classes `"coords"` with [`latlon<-()`](https://mark-eis.github.io/Waypoint/reference/latlon.html).

- Ensures values within objects of classes `"coords"` or `"waypoints"` are plausible geographic locations with `validate()`.

- Provides S3 `print()` methods for objects of classes `"coords"` or `"waypoints"` for neatly formatted printing.

*Waypoint* uses high performance C++ code seamlessly integrated into R using [`Rcpp`](https://www.rcpp.org).
