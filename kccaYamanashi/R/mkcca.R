
################################################################################
# Implementation of kernel cca method from "Extraction of correlated gene
# clusters from multiple genomic data by generalized kernel canonical
# correlation analysis." This is very similar to the kcca method in kernlab,
# except it let's us work with more than two tables at a time.
################################################################################


#' @title Merge default options for kernel cca
#' @param opts A partially specified list of options to use in kcca. The current
#' possible options are
#'   $kernels A list of kernlab kernel class objects. Defaults to a list of RBF
#'    kernels with bandwidth 1.
#'   $lambdas The regularization parameters for each table. Defaults to 1.
#' @param p The number of tables in the kernel CCA.
#' @return The original opts, but with unspecified options filled in with
#' defaults.
#' @importFrom kernlab rbfdot
#' @export
merge_kernel_opts <- function(opts = list(), p) {
  default_opts <- list()
  default_opts$kernels <- replicate(p, rbfdot(1), simplify = F)
  default_opts$lambdas <- rep(1, p)
  modifyList(default_opts, opts)
}

#' @title Get the kernel matrices associated with each table.
#' @param x_list A list whose elements are matrices with samples on rows and
#' features on columns.
#' @param opts A list of options specifying how to perform the kcca. See
#' merge_kernel_opts.
#' @return A list whose i^th element is the n x n kernel matrix associated with
#' the i^th data set.
#' @importFrom kernlab kernelMatrix
#' @export
get_k_matrices <- function(x_list, opts) {
  p <- length(x_list)
  k_list <- vector(length = p, mode = "list")
  for(i in seq_along(k_list)) {
    k_list[[i]] <- kernelMatrix(opts$kernels[[i]], x_list[[i]])
  }
  k_list
}

#' @title Get the left hand side matrix in equation (16) of Yamanashi's kcca
#' method
#' @param k_list A list whose i^th element is the n x n kernel matrix associated
#' with the i^th data set.
#' @param opts A list of options specifying how to perform the kcca. See
#' merge_kernel_opts.
#' @return The the left hand side matrix in equation (16) of Yamanashi's kcca
#' method.
#' @references "Extracting of correlated gene clusters from multiple genomic
#' data by generalized kernel canonical correlation analysis."
#' @export
mkcca_eigen_a <- function(k_list, opts) {
  n <- nrow(k_list[[1]])

  # A in eigenproblem Av = mu * B * v
  A <- matrix(0, n * p, n * p)
  for(i in seq_len(p)) {
    cur_i_ix <- (n * (i - 1) + 1) : (n * i)
    for(j in seq_len(i - 1)) {
      cur_j_ix <- (n * (j - 1) + 1) : (n * j)
      A[cur_i_ix, cur_j_ix] <- k_list[[i]] * k_list[[j]]
    }
  }
  A + t(A)
}

#' @title Get the right hand side matrix in equation (16) of
#' Yamanashi's kcca method
#' @param k_list A list whose i^th element is the n x n kernel matrix associated
#' with the i^th data set.
#' @param opts A list of options specifying how to perform the kcca. See
#' merge_kernel_opts.
#' @return The the right hand side matrix in equation (16) of Yamanashi's kcca
#' method.
#' @references "Extracting of correlated gene clusters from multiple genomic
#' data by generalized kernel canonical correlation analysis."
#' @export
mkcca_eigen_b <- function(k_list, opts) {
  n <- nrow(k_list[[1]])

  # B in eigenproblem Av = mu * B * v
  B <- matrix(0, n * p, n * p)
  for(i in seq_len(p)) {
    cur_ix <- (n * (i - 1) + 1) : (n * i)
    B[cur_ix, cur_ix] <- crossprod(k_list[[i]] + opts$lambdas[i] * diag(n))
  }
  B
}

#' @title Perform Yamanashi's multiple kernel canonical correlation analysis
#' @param k_list A list whose i^th element is the n x n kernel matrix associated
#' with the i^th data set
#' @param opts A list of options specifying how to perform the kcca. See
#' merge_kernel_opts.
#' @return The eigenvalues and eigenvectors for the eigenproblem (16) in the
#' referenced paper.
#' @references "Extracting of correlated gene clusters from multiple genomic
#' data by generalized kernel canonical correlation analysis."
#' @export
mkcca <- function(x_list, opts) {
  opts <- merge_kernel_opts(opts, p)
  k_list <- get_k_matrices(x_list, opts)
  A <- mkcca_eigen_a(k_list, opts)
  B <- mkcca_eigen_b(k_list, opts)
  eigen(solve(B) %*% A)
}
