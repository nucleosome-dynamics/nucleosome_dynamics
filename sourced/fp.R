#!/usr/bin/Rscript

# Some functional programming definitions

myFilter <- function (x, f, ...)
    # Vectorized version of Reduced made to work nicely with apply-family
    # functions
    x[f(x, ...)]

partial <- function (f, ...)
{   # Partial function application.
    # partial(f, x)(y) should be equivalent to f(y, x)
    capture <- list(...)
    function(x) do.call(f, c(list(x), capture))
}

compose <- function (...)
{   # Function composition.
    # compose(f, g, h)(x) sholud be equivalent to f(g(h(x)))
    comp2 <- function(f, g) {
        force(f)
        force(g)
        function(x) f(g(x))
    }
    Reduce(comp2, list(...))
}

flip2args <- function (f)
    # Flip the order of the arguments of a function that has two arguments
    function(x, y) f(y, x)

application <- function (f, ...)
    # Function application.
    # Useful to apply a list of functions on one same argument
    f(...)

procFun <- function (f, g)
{
    force(g)
    f(g)
}
