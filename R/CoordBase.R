# Waypoint R Package
# Mark Eisler Feb 2025
# For conversion and validation of geographic coordinates
#
# Requires R version 4.4.2 (2024-10-31) -- "Pile of Leaves" or later
#
# CoordBase.R

# ========================================
#  Validate Coordinates and Waypoints
#  S3generic validate()
#
#' @export

validate <- function(object, ...) 
    UseMethod("validate")
    

# ========================================
#  Analyse Coordinates and Waypoints
#  S3generic checkvalid()
#
#' @export

checkvalid <- function(object, ...) 
    UseMethod("checkvalid")


# ========================================
#  Analyse Coordinates
#  S3method checkvalid.coords()
#'
#' @export

checkvalid.coords <- function(coord, show_n = 20) {
    if (!inherits(coord, "coords"))
        stop("Argument `coords` must have class `\"coords\"`\n", call. = FALSE)
    invalid <- !attr(coord, "valid")
    n_invalid = sum(invalid, na.rm = TRUE)
    last_invalid = which(invalid)[show_n]
    if (n_invalid > show_n) {
        warning(
            n_invalid, " invalid coords, showing first ", show_n,
            "\n\t(use arg `show_n` to see more) ",
            call.= FALSE
        )
        tmp <- lapply(attributes(coord), \(x) x[seq_len(min(length(x), last_invalid))])
        coord <- coord[seq_len(last_invalid)]
        attributes(coord) <- tmp
        return(checkvalid(coord, show_n))
    }
    fmt <- attr(coord, "fmt");
    if (n_invalid) {
        invalids <- coord[invalid]
        suppressWarnings(coords(invalids) <- fmt)
        acl <- attr(coord, "latlon")
        if (!is.null(acl)) {
            if (length(acl) > 1)
                suppressWarnings(latlon(invalids) <- acl[invalid])
            else
                suppressWarnings(latlon(invalids) <- acl[1])
        }
    } else
        invalids <- NA_integer_
    list(
        allvalid = all(attr(coord, "valid")),
        n_invalid = n_invalid,
        invalids = invalids,
        which_invalid = which(invalid) %L% NA_integer_
    )
}

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
#' @param x,y atomic vector arguments or other objects for which `length()` is defined.
#'
#' @return \code{x}, or if \code{length(x)} is zero, \code{y}.
#'
#' @keywords univar arith
#' @export
#' @examples
#' c1 <- letters[1:4]
#' c2 <- character(0)
#' n1 <- 1:3
#' n2 <- numeric(0)
#' 
#' c1 %L% n1
#' c2 %L% n1
#' 
#' n1 %L% c1
#' n2 %L% c1
#'
#' rm(c1, c2, n1, n2)

`%L%` <- function (x, y) 
if (length(x)) x else y

#' Work in progress
#'

which_ll <- function(candidate = c("latitude", "longitude"))
{
    candidate <- tolower(candidate)
    match.arg(candidate)
}

