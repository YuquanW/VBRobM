% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bisquare_loss.R
\name{bisquare}
\alias{bisquare}
\alias{psi_bisquare}
\alias{dpsi_bisquare}
\alias{weight_bisquare}
\title{Bisquare Functions}
\usage{
psi_bisquare(x, k)

dpsi_bisquare(x, k)

weight_bisquare(x, k)
}
\arguments{
\item{x}{a vector of residuals to be robustfied.}

\item{k}{the tuning parameter in bisquare function.}
}
\value{
\code{psi_bisquare} returns \eqn{\psi(x)=\rho_{bisquare}^{\prime}(x)/2}, \code{dpsi_bisquare} returns \eqn{\psi^{\prime}(x)} and \code{weight_bisquare} returns \eqn{\psi(x)/x}.
}
\description{
Bisquare Functions
}
