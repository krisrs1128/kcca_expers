
################################################################################
# Implementation of multiple / integrated kernel cca method from "Extraction of
# correlated gene clusters from multiple genomic data by generalized kernel
# canonical correlation analysis."
################################################################################

#' @title Get the left hand side matrix in equations (16 / 18) of Yamanashi's
#' kcca method
#' @param k_list A list whose i^th element is a list of  n x n kernel matrices
#' associated with the i^th data set. For kcca, only the first kernel matrix
#' for each data set will be used.
#' @param opts A list of options specifying how to perform the kcca. See
#' merge_kernel_opts.
#' @return The the left hand side matrix in equation (16) of Yamanashi's kcca
#' method.
#' @references "Extracting of correlated gene clusters from multiple genomic
#' data by generalized kernel canonical correlation analysis."
#' @export
kcca_eigen_a <- function(k_list, opts) {
  n <- nrow(k_list[[1]])
  p <- length(k_list)

  # A in eigenproblem Av = mu * Bv
  A <- matrix(0, n * p, n * p)
  for(i in seq_len(p)) {
    cur_i_ix <- (n * (i - 1) + 1) : (n * i)
    for(j in seq_len(p)) {
      cur_j_ix <- (n * (j - 1) + 1) : (n * j)
      A[cur_i_ix, cur_j_ix] <- k_list[[i]] %*% k_list[[j]]
    }
  }
  A
}

#' @title Get the right hand side matrix in equations (16 / 18) of
#' Yamanashi's kcca method
#' @param k_list A list whose i^th element is a list of  n x n kernel matrices
#' associated with the i^th data set. For kcca, only the first kernel matrix
#' for each data set will be used.
#' @param opts A list of options specifying how to perform the kcca. See
#' merge_kernel_opts.
#' @return The the right hand side matrix in equation (16) of Yamanashi's kcca
#' method.
#' @references "Extracting of correlated gene clusters from multiple genomic
#' data by generalized kernel canonical correlation analysis."
#' @export
kcca_eigen_b <- function(k_list, opts) {
  n <- nrow(k_list[[1]])
  p <- length(k_list)

  # B in eigenproblem Av = mu * B * v
  B <- matrix(0, n * p, n * p)
  for(i in seq_len(p)) {
    cur_ix <- (n * (i - 1) + 1) : (n * i)
    B[cur_ix, cur_ix] <- crossprod(k_list[[i]] + opts$lambdas[i] * diag(n))
  }
  B
}

#' @title Perform Yamanashi's multiple / integrated kernel canonical correlation
#' analysis
#' @description Both multiple and integrated kccas are handled simultaneously,
#' by assuming the input data is a list whose i^th element is a list of kernels
#' to be applied to the i^th data set. The sum of these matrices is used in the
#' final eigenproblem, as in equation (18) of the paper. However, we can have
#' multiple attributes for multiple data sets, and each of the i element of the
#' data list will be expanded into a separate kcca block, as in equation (16).
#' That is, to do integrated kcca, ensure the $kernels element of opts is a list
#' of length 2, both of which are lists of kernels to apply to each data set. On
#' the other hand, to do multiple kcca, ensure that the $kernels element of opts
#' is a list of length 1 lists, containing the single kernel to apply to that
#' data set.
#' @param x_list A list whose elements are matrices with samples on rows and
#' features on columns.
#' @param opts A list of options specifying how to perform the kcca. See
#' merge_kernel_opts.
#' @return The eigenvalues and eigenvectors for the eigenproblems (16 / 18) in
#' the referenced paper.
#' @references "Extracting of correlated gene clusters from multiple genomic
#' data by generalized kernel canonical correlation analysis."
#' @importFrom magrittr %>%
#' @importFrom geigen geigen
#' @examples
#' # integrated kcca
#' x_list <- replicate(2, matrix(rnorm(100), 10, 10), simplify = F)
#' rbf_kernels <- list(rbfdot(1), rbfdot(10))
#' opts <- list(kernels = replicate(2, rbf_kernels, simplify = F))
#' kcca(x_list, opts)
#'
#' # multiple kcca
#' x_list <- replicate(10, matrix(rnorm(100), 10, 10), simplify = F)
#' opts <- list(kernels = replicate(10, list(rbfdot(1)), simplify = F))
#' kcca(x_list, opts)
#' @export
kcca <- function(x_list, opts) {
  opts <- merge_kernel_opts(opts, length(x_list))
  k_list <- get_k_matrices(x_list, opts)
  k_list_sums <- lapply(k_list, function(k) { center_mat(Reduce("+", k)) })
  A <- kcca_eigen_a(k_list_sums, opts)
  B <- kcca_eigen_b(k_list_sums, opts)
  eigen_res <- geigen(A, B)
  postprocess_res(eigen_res, k_list_sums)
}

#' @title Extract KCCA Scores
#' @param eigen_res The result of a call to geigen
#' @param k_list A list whose i^th element is a list of  n x n kernel matrices
#' @return A list with the following elements,
#'   $scores: A list whose i^th element is a matrix of scores for the i^th table. \cr
#'   $values: The eigenvalues of the problem, in descending order. \cr
#'   $vectors: The eigenvectors of the KCCA problem (not directly interpretable) \cr
#' @export
postprocess_res <- function(eigen_res, k_list) {
  n <- nrow(k_list[[1]])
  p <- length(k_list)
  scores <- vector(length = p, mode = "list")
  eigen_res$values <- rev(eigen_res$values)
  eigen_res$vectors <- eigen_res$vectors[, (n * p):1]

  for(i in seq_len(p)) {
    cur_ix <- (n * (i - 1) + 1) : (n * i)
    scores[[i]] <-  k_list[[i]] %*% eigen_res$vectors[cur_ix, ]
  }
  list(scores = scores, values = eigen_res$values, vectors = eigen_res$vectors)
}
