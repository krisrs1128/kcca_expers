
################################################################################
# Implementation of multiple kernel cca method from "Extraction of correlated
# gene clusters from multiple genomic data by generalized kernel canonical
# correlation analysis."
################################################################################

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
#' @importFrom magrittr %>%
#' @export
mkcca <- function(x_list, opts) {
  opts <- merge_kernel_opts(opts, p)
  k_list <- get_k_matrices(x_list, opts) %>%
    lapply(function(x) x[[1]])
  A <- mkcca_eigen_a(k_list, opts)
  B <- mkcca_eigen_b(k_list, opts)
  eigen(solve(B) %*% A)
}
