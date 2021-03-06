% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/kcca.R
\name{kcca_eigen_a}
\alias{kcca_eigen_a}
\title{Get the left hand side matrix in equations (16 / 18) of Yamanashi's
kcca method}
\usage{
kcca_eigen_a(k_list, opts)
}
\arguments{
\item{k_list}{A list whose i^th element is a list of  n x n kernel matrices
associated with the i^th data set. For kcca, only the first kernel matrix
for each data set will be used.}

\item{opts}{A list of options specifying how to perform the kcca. See
merge_kernel_opts.}
}
\value{
The the left hand side matrix in equation (16) of Yamanashi's kcca
method.
}
\description{
Get the left hand side matrix in equations (16 / 18) of Yamanashi's
kcca method
}
\references{
"Extracting of correlated gene clusters from multiple genomic
data by generalized kernel canonical correlation analysis."
}

