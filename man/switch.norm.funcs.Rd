% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/switch.norm.funcs.R
\name{switch.norm.funcs}
\alias{switch.norm.funcs}
\title{Switch Normalization Functions}
\usage{
switch.norm.funcs(norm.list = c("NN", "MN", "QN", "VSN"), norm.funcs = NULL)
}
\arguments{
\item{norm.list}{Switch all the build-in normalization methods into function names, including "NN", "MN", "QN", and "VSN".}

\item{norm.funcs}{New functions that user can create by themselves.}
}
\value{
A list that transforms the normalization method name into function name.
}
\description{
Switch the normalization method name into function name, for running in the normalization procedure
}
\details{
switch.norm.funcs
}
\examples{
switch.norm.funcs(norm.list = c("NN", "MN", "QN"), norm.funcs = NULL)

}
\keyword{names}
\keyword{switch}
