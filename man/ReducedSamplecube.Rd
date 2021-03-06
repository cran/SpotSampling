% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReducedSamplecube.R
\name{ReducedSamplecube}
\alias{ReducedSamplecube}
\title{Cube method with reduction of the auxiliary variables matrix}
\usage{
ReducedSamplecube(X, pik, redux = TRUE, t)
}
\arguments{
\item{X}{a matrix of size (N x p) of auxiliary variables on which the sample must be balanced.}

\item{pik}{a vector of size N of inclusion probabilities.}

\item{redux}{a boolean value that specify if matrix \code{X} is reduced during the cube method. Default value is TRUE.}

\item{t}{the maximum number of constraints that can potentially be removed during the landing phase.}
}
\value{
the updated vector of \code{pik} that contains only 0s and 1s that indicates if a unit is selected or not at each wave.
}
\description{
Modified cube method.
This function reduces considerably the execution time when the matrix of auxiliary variables \code{X} contains lot of 0s.
It is based on the function \code{\link[sampling:samplecube]{samplecube}} from the package \code{sampling}.
}
\details{
In case where the number of auxiliary variables is great (i.e. p very large), even if we use the fast implementation proposed by
(Chauvet and Tille 2005), the problem is time consuming.
This function reduces considerably the execution time when the matrix of auxiliary variables \code{X} contains lot of 0s.
It considers a reduced matrix \code{X} by removing columns and rows that sum to 0 (see  \code{\link{ReducedMatrix}}).

Moreover, the landing by variable suppression is used.
\code{t} specifies the maximum number of constraints that can potentially be removed during the landing phase. This means that the first (N-T) constraints in \code{X} can be exactly satisfied.
}
\examples{
set.seed(1)
## Matrix of 8 auxilary variables and 10 units with lot of 0s ##
X <- matrix(c(0.6,0.0,0.0,0.0,
              0.1,0.0,0.1,0.0,
              0.3,0.0,0.0,0.3,
              0.0,0.3,0.0,0.3,
              0.0,0.6,0.0,0.0,
              0.0,0.1,0.1,0.0), ncol = 4, byrow = TRUE)

## Inclusion probabilities with 10 units ##
pik <- c(0.60,0.10,0.30,0.30,0.60,0.10)

## parameter t ##
t <- 2

## Cube method ##
s   <- ReducedSamplecube(X, pik, redux = TRUE, t)
s

}
\references{
Chauvet, G. and Tille, Y. (2006). A fast algorithm of balanced sampling. Computational Statistics, 21/1:53-62
}
\seealso{
\code{\link[sampling:samplecube]{samplecube}}, \code{\link[sampling:landingcube]{landingcube}}, \code{\link{ReducedMatrix}}.
}
\author{
Esther Eustache \email{esther.eustache@unine.ch}, Raphael Jauslin \email{raphael.jauslin@unine.ch}
}
