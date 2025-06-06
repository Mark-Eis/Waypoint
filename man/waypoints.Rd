% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoordBase.R, R/RcppExports.R
\name{waypoints}
\alias{waypoints}
\alias{as_waypoints}
\alias{as_waypoints.default}
\title{Geographic Waypoint Class}
\usage{
as_waypoints(object, ...)

\method{as_waypoints}{default}(object, ..., fmt = 1L)
}
\arguments{
\item{object}{a data frame with each row representing a waypoint, comprising at least two
\code{numeric} columns containing values of latitude and longitude, and optionally a
\code{character} column of waypoint names (see \emph{Details}).}

\item{...}{further arguments passed to or from other methods.}

\item{fmt}{\code{integer}, \code{1L}, \code{2L} or \code{3L}, specifying the required coordinate
format.}
}
\value{
An object of classes \code{"waypoints"} and \code{"data.frame"}, comprising the original data
frame argument \code{object}, with additional attributes: –

\item{\code{"class"}}{the \code{character} string "waypoints".}

\item{\code{"fmt"}}{an \code{integer} indicating the coordinate format.}

\item{\code{"namescol"}}{an \code{integer} indicating the position of a waypoint names column,
  if present.}

\item{\code{"llcols"}}{a length 2 \code{integer} vector indicating the positions of latitude and
  longitude columns.}

\item{\code{"validlat"} and \code{"validlon"}}{\code{logical} vectors indicating whether
  individual latitude and longitude values are valid geographic locations.}
}
\description{
\code{as_waypoints()} creates an object of class \code{"waypoints"}, a robust representation of a
series of geographic or GPS waypoints of paired latitude and longitude values.
}
\details{
Data frame argument \code{object} should have \code{numeric} vector latitude and longitude
columns with individual values having a decimal point after the number of whole degrees in the
case of \emph{decimal degrees}, after the number of whole minutes in the case of
\emph{degrees and minutes}, and after the number of whole seconds in the case of
\emph{degrees, minutes and seconds}. These should be the first two columns of the data frame, or
the second and third columns if the first column contains waypoints names (see below).
Alternative columns may be specified for the latitude and longitude by setting \code{"llcols"} as
a length 2 \code{integer} vector attribute of \code{object} indicating their positions in the
data frame.

The \code{fmt} argument should be \code{1L} to represent decimal degrees, \code{2L} for degrees
and minutes, and \code{3L} for degrees, minutes and seconds and is used to provide the
format of values in data frame argument \code{object} to be converted to class
\code{"waypoints"}.

If the waypoints have names, these should be included in a "Name" column of data frame argument
\code{object}, by default immediately before (on the left-hand side of) the latitude and
longitude columns. An alternative column for waypoint names may be specified by setting an
\code{integer} \code{\link[base:attributes]{attribute}} named "namescol" indicating its position
in \code{object}. If neither a "Name" column nor a \code{"namescol"} attribute is present in
\code{object}, \code{\link[base:row.names]{row.names}} are used for waypoint names.

The latitude and longitude values of a newly created \code{"waypoints"} object are checked to
ensure they are valid geographic locations as described under \code{\link{validate}()}.
}
\examples{
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

## Show as an ordinary R data frame
as.data.frame(wp1)

###

## Dataframe representing unnamed latitude and longitude
## values in decimal degrees
wp2 <- data.frame(
    lat = c(51.507765, 49.54621, 48.107232, 38.889494, 0, -37.11174, -53.104781, -25.240156),
    lon = c(-0.127924, 18.398562, -122.778671, -77.035242, 0, -12.28863, 73.517283, -57.519227)
)

## Create unnamed "waypoints" object in decimal degrees (default fmt = 1)
as_waypoints(wp2)

## Add waypoint names as row.names
row.names(wp2) <-
    c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
      "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport")

wp2

## Show as an ordinary R data frame
as.data.frame(wp2)

rm(wp1, wp2)


}
\seealso{
\code{\link[base:attr]{attr}()}, \code{\link[base:data.frame]{data.frame}()}, and
  \code{\link{validate}()}.

Other coordsandway: 
\code{\link{convert}()},
\code{\link{coords}}
}
\concept{coordsandway}
