# Operator Providing Alternative to Zero-Length Object

Infix function implementing provision of an alternative if an object has
zero length.

## Usage

``` r
x %L% y
```

## Arguments

- x, y:

  atomic vector arguments or other objects for which
  [`length()`](https://rdrr.io/r/base/length.html) is defined.

## Value

`x`, or if `length(x)` is zero, `y`.

## Details

The infix function `%L%` may be useful in implementing
`if (length(x)) x else y` and was inspired by the null coalescing
operator [`%||%`](https://rdrr.io/r/base/Control.html).

## See also

[`%||%`](https://rdrr.io/r/base/Control.html).

## Examples

``` r
c4 <- letters[1:4]
c0 <- character(0)
n3 <- 1:3
n0 <- numeric(0)

c4 %L% n3
#> [1] "a" "b" "c" "d"
c0 %L% n3
#> [1] 1 2 3

n3 %L% c4
#> [1] 1 2 3
n0 %L% c4
#> [1] "a" "b" "c" "d"

rm(c4, c0, n3, n0)
```
