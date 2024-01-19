#' Bisquare Functions
#'
#' @param x a vector of residuals to be robustfied.
#' @param k the tuning parameter in bisquare function.
#'
#' @return `psi_bisquare` returns \eqn{\psi(x)=\rho_{bisquare}^{\prime}(x)/2}, `dpsi_bisquare` returns \eqn{\psi^{\prime}(x)} and `weight_bisquare` returns \eqn{\psi(x)/x}.
#'
#' @name bisquare
NULL

#' @rdname bisquare
#' @export
psi_bisquare <- function(x, k) {
  #return(x*(abs(x) <= k) + sign(x)*k*(abs(x) > k))
  return(x*(1 - (x/k)^2)^2*(abs(x) <= k))
}

#' @rdname bisquare
#' @export
dpsi_bisquare <- function(x, k) {
  #return(1*(abs(x) <= k) + 0*(abs(x) > k))
  return(((1 - (x/k)^2)^2 - 4*x^2/k^2*(1 - (x/k)^2))*(abs(x) <= k))
}

#' @rdname bisquare
#' @export
weight_bisquare <- function(x, k) {
  #return(pmin(1, k/abs(x)))
  return((1 - (x/k)^2)^2*(abs(x) <= k))
}


