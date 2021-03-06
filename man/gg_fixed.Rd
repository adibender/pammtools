% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convenience-plots.R
\name{gg_fixed}
\alias{gg_fixed}
\title{Forrest plot of fixed coefficients}
\usage{
gg_fixed(x, intercept = FALSE, ...)
}
\arguments{
\item{x}{A model object.}

\item{intercept}{Logical, indicating whether intercept term should be included.
Defaults to \code{FALSE}.}

\item{...}{Currently not used.}
}
\description{
Given a model object, returns a data frame with columns \code{variable},
\code{coef} (coefficient), \code{ci_lower} (lower 95\\% CI) and
\code{ci_upper} (upper 95\\% CI).
}
\examples{
g <- mgcv::gam(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species,
 data=iris)
gg_fixed(g, intercept=TRUE)
gg_fixed(g)
}
\seealso{
\code{\link{tidy_fixed}}
}
