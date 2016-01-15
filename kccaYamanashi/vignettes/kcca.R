
################################################################################
# R script accompanying kcca.rnw
################################################################################

## ---- libraries ----
library("ggplot2")
library("kernlab")
library("kccaYamanashi")
library("reshape2")

## ---- simulate-data ----
n <- 40
C <- sapply(1:2, function(x) sample(0:1, n, replace = T))
theta <- C + matrix(runif(n * 2, 0, 1), n, 2)
theta <- theta / rowSums(theta)
theta <- rbeta(n, 2, 2)
theta <- cbind(theta, 1 - theta)

p <- 1000
x <- seq(-4, 4, length.out = p)
g <- list(dnorm(x, -1), dnorm(x, 1))

f <- t(apply(theta, 1, function(z) z[1] * g[[1]] + z[2] * g[[2]]))

## ---- vis-simulated-data ----
mf <- melt(f)
mf$x <- x[mf$Var2]
mf$theta <- theta[, 1]
ggplot(mf) +
  geom_line(aes(x = x, y = value, group = Var1, col = theta))

## ---- run-kcca ----
sigma <- 3
kernels <- list(list(rbfdot(sigma)), list(vanilladot()))
kcca_res <- kccaYamanashi::kcca(list(f, theta), list(kernels = kernels))

## ---- eivenvals ----
ggplot(data.frame(ix = 1:10, lambda = kcca_res$values[1:10])) +
  geom_bar(aes(x = ix, y = lambda), stat = "identity") +
  ggtitle("top eigenvalues")

## ---- kcca-scores-2 ----
u <- kcca_res$scores[[2]][, 1:2]
ggplot(data.frame(theta = theta[, 1], u = u)) +
  geom_point(aes(x = u.1, y = u.2, col = theta)) +
  ggtitle("kcca scores")

## ---- correlation-scores ----
ggplot(data.frame(theta = theta[, 1], ux = u[1:n], uy = u[(n + 1):(2 * n)])) +
  geom_point(aes(x = ux, y = uy)) +
  ggtitle("Correlation between two table scores")

## ---- top-eigen ----
v <- kcca_res$vectors[, 1]
ggplot(data.frame(ix = 1:length(v), v, theta = theta[, 1])) +
  geom_point(aes(x = ix, y = v, col = theta)) +
  ggtitle("top eigenvector")

## ---- second-eigen ----
v <- kcca_res$vectors[, 2]
ggplot(data.frame(ix = 1:length(v), v, theta = theta[, 1])) +
  geom_point(aes(x = ix, y = v, col = theta)) +
  ggtitle("second eigenvector")

## ---- pca-approach ----
pca_res <- svd(f)
pca_res <- pca_res$u %*% diag(pca_res$d)
ggplot(data.frame(pca_res = pca_res, theta = theta[, 1])) +
  geom_point(aes(x = pca_res.1, y = pca_res.2, col = theta))
