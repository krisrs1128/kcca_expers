
################################################################################
# Utilities used in both multiple and integrated kernel cca approaches.
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
merge_kernel_opts <- function(opts = list(), p = 1) {
  default_opts <- list()
  default_opts$kernels <- replicate(p, list(rbfdot(1)), simplify = F)
  default_opts$lambdas <- rep(1, p)
  modifyList(opts, default_opts)
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
    cur_kernels <- vector(length = length(opts$kernels[[i]]), mode = "list")
    for(j in seq_along(opts$kernels[[i]])) {
      cur_kernels[[j]] <- kernelMatrix(opts$kernels[[i]][[j]], x_list[[i]])
    }
    k_list[[i]] <- cur_kernels
  }
  k_list
}

