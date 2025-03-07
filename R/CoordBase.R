## Waypoint R Package
## Mark Eisler (c) March 2025
## For conversion and validation of geographic coordinates
##
## Requires R version 4.4.2 (2024-10-31) -- "Pile of Leaves" or later
##
## CoordBase.R

## __________________________________________________
## Set R vector object class to coords and return,
## or convert format of R coords object and return
#' @title Geographic or GPS Coordinate Class
#' 
#' @name Coords
#' 
#' @description
#' \code{as_coords()} creates an object of class \code{"coords"}, a robust representation of a
#' series of geographic or GPS coordinate values.
#' 
#' \code{convert()} converts the format of existing objects of class \code{"coords"} between (i)
#' decimal degrees, (ii) degrees and minutes, and (iii) degrees, minutes and seconds.
#'
#' @details
#' Individual values provided in the \code{numeric} vector argument \code{nv} should have a decimal
#' point after the number of whole degrees in the case of \emph{decimal degrees}, after the number
#' of whole minutes in the case of \emph{degrees and minutes}, and after the number of whole
#' seconds in the case of \emph{degrees, minutes and seconds}.
#'
#' The \code{fmt} argument should be \code{1L} to represent decimal degrees, \code{2L} for degrees
#' and minutes, and \code{3L} for degrees, minutes and seconds and is used to provide both the
#' format of values in \code{numeric} vector argument \code{nv} to be converted into a
#' \code{"coords"} object and the desired format if a \code{"coords"} object is to be converted to
#' a new format. Note that on conversion of a \code{"coords"} object, the original \code{numeric}
#' vector argument \code{nv} is modified such that the values are as described in the previous
#' paragraph, and may be inspected using standard R code, see examples.
#'
#' The values of a newly created \code{"coords"} object are checked to ensure they are valid
#' geographic locations as described under \code{\link[=validate]{validate}()}. Likewise, a
#' check is made to ensure that an existing \code{"coords"} object to be converted to a new format
#' has already been validated; if not, it is re-validated. 
#'
#' @family coords_waypoints
#' @seealso
#' \code{\link[base:attr]{attr}()}, \code{\link[base:attributes]{attributes}},
#'   \code{\link[=latlon]{latlon}()}, \code{\link[base:numeric]{numeric}()} and
#'   \code{\link[=validate]{validate}()}.
#'
#' @param object a \code{numeric} vector of coordinate values, optionally named, or an object of
#'   class \code{"waypoints"}.
#'
#' @param \dots further arguments passed to or from other methods.
#'
#' @param x object of class \code{"coords"} created by function
#'   \code{\link[=as_coords]{as_coords}()}.
#'
#' @param fmt \code{integer}, 1L, 2L or 3L, indicating the current or desired coordinate format.
#'
#' @param usenames \code{logical}, whether or not to include names in formatted output.
#'
#' @param latlon \code{logical}, indicating whether the \code{as_coords()} S3 method for class
#'   \code{"waypoints"} extracts the latitude component of argument \code{object} (if \code{TRUE}),
#'   or the longitude (if \code{FALSE}).
#'
#' @return
#' An object of class \code{"coords"}, comprising the original a \code{numeric} vector argument
#' \code{nv} with values possibly converted as appropriate and additional attributes: –
#' \item{\code{"fmt"}}{the coordinate format.}
#' \item{\code{"valid"}}{a \code{logical} vector indicating whether individual coordinate values
#'   are valid geographic locations.}
#'
#' @examples
#' ## Numeric vector representing degrees and minutes
#' dm <- c(5130.4659, 4932.7726, 4806.4339, 3853.3696, 0.0000, -3706.7044, -5306.2869, -2514.4093,
#'         -007.6754, 1823.9137, -12246.7203, -7702.1145, 0.0000, -1217.3178, 7331.0370, -5731.1536)
#'
#' ## Create a "coords" object of degrees and minutes (fmt = 2)
#' as_coords(dm, fmt = 2)
#'
#' ## Name the "coords" object
#' names(dm) <- rep(c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
#'                    "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"), 2)
#' dm
#'
#' ## Convert to degrees, minutes and seconds (fmt = 3)
#' convert(dm, 3)
#'
#' ## Convert to decimal degrees (fmt = 1)
#' convert(dm, 1)
#'
#' ## Decimal degrees as an ordinary R numeric vector
#' as.numeric(dm)
#'
#' ## Convert to degrees and minutes, then format as a
#' ## fixed-width character vector without names...
#' convert(dm, 3) |> format(usenames = FALSE)
#'
#' ## ...or with them
#' format(dm)
#'
#' rm(dm)
#'

## ========================================
##  Create Coordinates
##  S3generic as_coords(object, ...)
#'
#' @export

as_coords <- function(object, ...) 
    UseMethod("as_coords")

## __________________________________________________
## Add "waypoints" to R data.frame object class and validate,
## or convert format of R waypoints object and return
#' @title Geographic or GPS Waypoint Class
#' 
#' @name Waypoints
#' 
#' @description
#' \code{as_waypoints()} creates an object of class \code{"waypoints"}, a robust representation of a
#' series of geographic or GPS waypoints of paired latitude and longitude values.
#' 
#' \code{convert()} converts the format of existing objects of class \code{"waypoints"} between (i)
#' decimal degrees, (ii) degrees and minutes, and (iii) degrees, minutes and seconds.
#'
#' @details
#' By default, the names of the waypoints should be included in a "Name" column of data frame
#' argument \code{df}, and the latitude and longitude in the two columns immediately on the right
#' hand side of "Name". An alternative column for waypoint names may be specified by setting an
#' \code{integer} attribute, \code{"namescol"} indicating its position in \code{df}, while setting
#' this attribute to \code{NA} supresses printing of waypoint names. If \code{df} has neither a
#' "Name" column nor a \code{"namescol"} attribute, the \code{"row.names"} attribute is used for
#' waypoint names if present in \code{df}. Similarly, alternative columns for the latitude and
#' longitude may be specified by setting \code{"llcols"} as a length 2 \code{integer} vector
#' attribute indicating their positions in \code{df}.
#'
#' Individual values provided in the \code{numeric} vector latitude and longitude columns of data
#' frame argument \code{df} should have a decimal point after the number of whole degrees in the
#' case of \emph{decimal degrees}, after the number of whole minutes in the case of
#' \emph{degrees and minutes}, and after the number of whole seconds in the case of
#' \emph{degrees, minutes and seconds}.
#'
#' The \code{fmt} argument should be \code{1L} to represent decimal degrees, \code{2L} for degrees
#' and minutes, and \code{3L} for degrees, minutes and seconds and is used to provide both the
#' format of values in data frame argument \code{df} to be converted into a \code{"waypoints"}
#' object and the desired format if a \code{"waypoints"} object is to be converted to a new format.
#' Note that on conversion of a \code{"waypoints"} object, the original data frame argument
#' \code{df} is modified such that the latitude and longitude values are as described in the
#' previous paragraph, and may be inspected using standard R code, see examples.
#'
#' The latitude and longitude values of a newly created \code{"waypoints"} object are checked to
#' ensure they are valid geographic locations as described under
#' \code{\link[=validate]{validate}()}. Likewise, a check is made to ensure that an existing
#' \code{"waypoints"} object to be converted to a new format has already been validated; if not, it
#' is re-validated. 
#'
#' @family coords_waypoints
#' @seealso
#' \code{\link[base:attr]{attr}()}, \code{\link[base:attributes]{attributes}},
#'   \code{\link[base:data.frame]{data.frame}()}, and \code{\link[=validate]{validate}()}.
#'
#' @param object a data frame with each row representing a waypoint, comprising at least two
#'   \code{numeric} columns containing values of latitude and longitude, and optionally a
#'   \code{character} column of waypoint names (see \emph{Details}). 
#'
#' @param \dots further arguments passed to or from other methods.
#'
#' @param x an object of class \code{"waypoints"} created by function
#' \code{\link[=as_waypoints]{as_waypoints}()}.
#'
#' @param fmt an \code{integer} of value 1L, 2L or 3L, indicating the current or desired coordinate
#'   format (see \emph{Details}).
#'
#' @param usenames \code{logical}, whether or not to include waypoint names in formatted output.
#'
#' @return
#' An object of classes \code{"waypoints"} and \code{"data.frame"}, comprising the original data
#' frame argument \code{df}, with latitude and longitude values possibly converted as appropriate
#' and additional attributes: –
#' \item{\code{"fmt"}}{the coordinate format.}
#' \item{\code{"namescol"}}{the position of waypoint names, if present within \code{df}.}
#' \item{\code{"llcols"}}{the position of latitude and longitude columns within \code{df}.}
#' \item{\code{"validlat"} and \code{"validlon"}}{\code{logical} vectors indicating whether
#'   individual latitude and longitude values are valid geographic locations.}
#'
#' @examples
#' ## Dataframe representing waypoint names, and latitude and longitude values
#' ## of degrees, minutes and seconds
#' wp1 <- data.frame(
#'     name = c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
#'              "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"),
#'     lat = c(513027.95, 493246.36, 480626.04, 385322.18, 0, -370642.26, -530617.21, -251424.56),
#'     lon = c(-00740.53, 182354.82, -1224643.22, -770206.87, 0, -121719.07, 733102.22, -573109.21)
#' )
#'
#' ## Create "waypoints" object of degrees, minutes and seconds (fmt = 3)
#' as_waypoints(wp1, fmt = 3)
#'
#' ## Convert to degrees and minutes (fmt = 2)
#' convert(wp1, 2)
#'
#' ## Convert to decimal degrees (fmt = 1)
#' convert(wp1, 1)
#'
#' ###
#' ## Dataframe representing unnamed latitude and longitude
#' ## values in decimal degrees
#' wp2 <- data.frame(
#'     lat = c(51.507765, 49.54621, 48.107232, 38.889494, 0, -37.11174, -53.104781, -25.240156),
#'     lon = c(-0.127924, 18.398562, -122.778671, -77.035242, 0, -12.28863, 73.517283, -57.519227)
#' )
#'
#' ## Create "waypoints" object of decimal degrees (default fmt = 1)
#' as_waypoints(wp2, fmt = 1)
#'
#' ## Add row.names
#' row.names(wp2) <-
#'     c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
#'       "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport")
#' wp2
#'
#' ## Convert to degrees and minutes (fmt = 2)
#' convert(wp2, 2)
#'
#' ## Convert to degrees, minutes and seconds (fmt = 3)
#' convert(wp2, 3)
#'
#' ## Degrees, minutes and seconds values as an ordinary R data frame
#' as.data.frame(wp2)
#'
#' ## Convert to decimal degrees, then format as a
#' ## fixed-width character vector without names...
#' convert(wp2, 3) |> format(usenames = FALSE)
#'
#' ## ...or with them
#' format(wp2)
#'
#' rm(wp1, wp2)
#'
## ========================================
##  Create Waypoints
##  S3generic as_waypoints(object, ...)
#'
#' @export

as_waypoints <- function(object, ...) 
    UseMethod("as_waypoints")


## ========================================
##  Convert Coordinates and Waypoints
##  S3generic convert(x, ...)
#'
#' @rdname Coords
#' @export

convert <- function(x, ...) 
    UseMethod("convert")


## ========================================
##  Print Coordinates
##  S3method print.coords(x, ...)
#'
#' @rdname Coords
#' @export

print.coords <- function (x, ...) {
    writeLines(format(x, ...))
    invisible(x)
}


## ========================================
##  Print Waypoints
##  S3method print.waypoints(x, ...)
#'
#' @rdname Waypoints
#' @export

print.waypoints <- function (x, ...) {
    fmtx <- format(x, ...)
    writeLines(ll_headers(fmtx, attr(x, "fmt")))
    writeLines(fmtx)
    invisible(x)
}


## ========================================
##  Validate Coordinates and Waypoints
##  S3generic validate(x, ...)
##
#' @export

validate <- function(x, ...) 
    UseMethod("validate")
    
## ========================================
#' @title
#' Review Coordinates and Waypoints Validity
#'
#' @description
#' \code{review()} review validity of elements of \code{"coords"} and \code{"waypoints"} objects.
#'
#' @details
#' \code{review()} reveals elements of \code{"coords"} and  \code{"waypoints"} objects that do not
#' conform to the criteria checked by \code{\link[=validate]{validate}()}, i.e. are not valid
#' geographic locations.
#'
#' @family validate
#' @seealso
#' \code{"\link[=as_coords]{coords}"} and \code{"\link[=as_waypoints]{waypoints}"}.
#'
#' @param x object of class \code{"coords"} or \code{"waypoints"}.
#'
#' @param show_n \code{integer}, the maximum number of invalid elements of argument \code{x} to
#' include in the output; default \code{20L}.
#'
#' @param \dots further arguments passed to or from other methods.
#'
#' @return
#' The \code{review()} method for class \code{"coords"} returns a \code{\link[base:list]{list}} comprising the
#' following elements: -
#'
#' \item{allvalid}{\code{logical}, whether or not all the elements of argument \code{x} are valid.}
#' \item{n_invalid}{\code{integer}, the number of invalid elements in argument \code{x}, if any.}
#' \item{invalids}{\code{numeric} vector including invalid elements of argument \code{x}, if any.}
#' \item{which_invalid}{\code{integer} vector specifying which elements of argument \code{x} are
#' 	 invalid, if any.}
#'
#' The method for class \code{"waypoints"} returns a list of two sub-lists, each sub-list with
#' elements as described for the method for class \code{"coords"}, one each for latitude and
#' longitude.
#'
#' @export
#'
#' @examples
#' ## Continuing example from `validate()`, named numeric vector representing degrees and minutes,
#' ## the erroneous first value having more than 60 minutes
#' \dontshow{
#'    oldopt <- options(warn = -1)
#'    dm <-
#'        c(5160.4659, 4932.7726, 4806.4339, 3853.3696, 0.0000, -3706.7044, -5306.2869, -2514.4093,
#'		   -007.6754, 1823.9137, -12246.7203, -7702.1145, 0.0000, -1217.3178, 7331.0370, -5731.1536)
#'
#'    names(dm) <- 
#'        rep(c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
#'              "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"), 2)
#' }
#'
#' ## Create "coords" object of degrees and minutes
#' as_coords(dm, fmt = 2)
#'
#' review(dm)
#'
#' ###
#' ## Continuing example from `validate()`, data frame representing waypoint names and latitude
#' ## and longitude values in decimal degrees, the erroneous penultimate latitude having more than
#' ## 90 degrees absolute value
#' \dontshow{
#' wp1 <- data.frame(
#'     name = c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
#'              "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"),
#'     lat = c(51.507765, 49.54621, 48.107232, 38.889494, 0, -37.11174, -93.104781, -25.240156),
#'     lon = c(-0.127924, 18.398562, -122.778671, -77.035242, 0, -12.28863, 73.517283, -57.519227)
#' )
#' }
#'
#' ## Create "waypoints" object of decimal degrees (fmt = 1)
#' as_waypoints(wp1, fmt = 1)
#'
#' review(wp1)
#'
#' rm(dm, wp1)
#' 
#' \dontshow{
#'    options(oldopt)
#' }

## ========================================
##  Review Coordinates and Waypoints
##  S3generic review(x, ...)
##
review <- function(x, ...) 
    UseMethod("review")

## ========================================
##  Review Coordinates
##  S3method review.coords(x, ..., show_n = 20L)
#'
#' @rdname review
#' @export

review.coords <- function(x, ..., show_n = 20L)
{
    if (!inherits(x, "coords"))
        stop("Argument `coords` must have class `\"coords\"`\n", call. = FALSE)
    invalid <- !attr(x, "valid")
    n_invalid = sum(invalid, na.rm = TRUE)
    last_invalid = which(invalid)[show_n]
    if (n_invalid > show_n) {
        warning(
            n_invalid, " invalid coords, showing first ", show_n,
            "\n\t(use arg `show_n` to see more) ",
            call.= FALSE
        )
        tmp <- lapply(attributes(x), \(x) x[seq_len(min(length(x), last_invalid))])
        x <- x[seq_len(last_invalid)]
        attributes(x) <- tmp
        return(review(x, show_n))
    }
    fmt <- attr(x, "fmt");
    if (n_invalid) {
        invalids <- x[invalid]
        suppressWarnings(as_coords(invalids, fmt = fmt))
        acl <- attr(x, "latlon")
        if (!is.null(acl)) {
            if (length(acl) > 1)
                suppressWarnings(latlon(invalids) <- acl[invalid])
            else
                suppressWarnings(latlon(invalids) <- acl[1])
        }
    } else
        invalids <- NA_integer_
    list(
        allvalid = all(attr(x, "valid")),
        n_invalid = n_invalid,
        invalids = invalids,
        which_invalid = which(invalid) %L% NA_integer_
    )
}


## ========================================
##  Review Waypoints
##  S3method review.waypoints(x, ..., show_n = 20L)
#'
#' @rdname review
#' @export

review.waypoints <- function(x, ..., show_n = 20L)
	lapply(c(TRUE, FALSE), \(y) review(as_coords(x, y), show_n)) |> setNames(c("Lat", "Lon"))


## ========================================
#' @title
#' Operator Providing Alternative to Zero-Length Object
#'
#' @name op-zero-length
#'
#' @description
#' Infix function implementing provision of an alternative if an object has zero length.
#'
#' @details
#' The infix function \code{\%L\%} may be useful in implementing \code{if (length(x)) x else y} and was inspired by
#' the null coalescing operator \code{\link[base:Control]{\%||\%}}.
#'
#' @family utils
#' @seealso \code{\link[base:Control]{\%||\%}}.
#'
#' @param x,y atomic vector arguments or other objects for which \code{length()} is defined.
#'
#' @return \code{x}, or if \code{length(x)} is zero, \code{y}.
#'
#' @keywords univar arith
#' @export
#' @examples
#' c4 <- letters[1:4]
#' c0 <- character(0)
#' n3 <- 1:3
#' n0 <- numeric(0)
#' 
#' c4 %L% n3
#' c0 %L% n3
#' 
#' n3 %L% c4
#' n0 %L% c4
#'
#' rm(c4, c0, n3, n0)

`%L%` <- function (x, y) 
if (length(x)) x else y

