# Waypoint R Package
# Mark Eisler Feb 2025
# For conversion and validation of geographic coordinates
#
# Requires R version 4.4.2 (2024-10-31) -- "Pile of Leaves" or later
#
# CoordBase.R

# ========================================
#  Validate Coordinates and Waypoints
#  S3generic validate(x, ...)
#
#' @export

validate <- function(x, ...) 
    UseMethod("validate")
    

# ========================================
#' @title
#' Review Coordinates and Waypoints
#'
#' @description
#' Review validity of elements of \code{"coords"} and  \code{"waypoints"} objects
#'
#' @details
#' \code{review()} reveals elements of \code{"coords"} and  \code{"waypoints"} objects that do not
#' conform to the criteria checked by \code{\link[=validate]{validate}()}, i.e. are not valid
#' geographic locations.
#'
#' @family validate
#' @seealso
#' \code{\link[=coords]{"coords"}} and \code{\link[=waypoints]{"waypoints"}}.
#'
#' @param x object of class \code{"coords"} or \code{"waypoints"}.
#'
#' @param show_n \code{integer}, the maximum number of invalid elements of argument \code{x} to
#' include in the output; default \code{20L}.
#'
#' @param \dots further arguments passed to or from other methods.
#'
#' @return
#' A \code{list} comprising the following components: -
#'
#' \item{allvalid}{\code{logical}, whether or not all the elements of argument \code{x} are valid.}
#' \item{n_invalid}{\code{integer}, the number of invalid elements in argument \code{x}, if any.}
#' \item{invalids}{\code{numeric} vector including invalid elements of argument \code{x}, if any.}
#' \item{which_invalid}{\code{integer} vector specifying which elements of argument \code{x} are
#' 	 invalid, if any.}
#'
#' @export
#'
#' @examples
#' ## Continuing example from `coords()`, named numeric vector representing degrees and minutes,
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
#' coords(dm) <- 2
#'
#' review(dm)
#'
#' ###
#' ## Continuing example from `waypoints()`, data frame representing waypoint names and latitude
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
#' ## Create "waypoints" object of decimal degrees
#' waypoints(wp1) <- 1
#'
#' review(wp1)
#'
#' rm(dm, wp1)
#' 
#' \dontshow{
#'    options(oldopt)
#' }

# ========================================
#  Review Coordinates and Waypoints
#  S3generic review(x, ...)
#
review <- function(x, ...) 
    UseMethod("review")


# ========================================
#  Review Coordinates and Waypoints
#  S3method review.coords(x, ...)
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
        suppressWarnings(coords(invalids) <- fmt)
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


# ========================================
#  Review Coordinates and Waypoints
#  S3method review.waypoints(x, ...)
#'
#' @rdname review
#' @export

review.waypoints <- function(x, ..., show_n = 20L)
	lapply(c(TRUE, FALSE), \(y) review(as_coord(x, y), show_n)) |> setNames(c("Lat", "Lon"))

# ========================================
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

