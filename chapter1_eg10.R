# Load libraries
library(glmnet)

# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

n <- 500
x <- array(dim = c(n, 6))
z <- array(dim = c(n, 2))
for (i in 1:2) {
  z[, i] <- rnorm(n)
}
y <- 3 * z[, 1] - 1.5 * z[, 2] + 2 * rnorm(n)
for (j in 1:3) {
  x[, j] <- z[, 1] + rnorm(n) / 5
}
for (j in 4:6) {
  x[, j] <- z[, 2] + rnorm(n) / 5
}
best_score <- Inf
for (alpha in seq(0, 1, 0.01)) {
  res <- cv.glmnet(x, y, alpha = alpha)
  lambda <- res$lambda.min
  min_cvm <- min(res$cvm)
  if (min_cvm < best_score) {
    alpha_min <- alpha
    lambda_min <- lambda
    best_score <- min_cvm
  }
}
alpha_min
lambda_min
glmnet(x, y, alpha = alpha_min, lambda = lambda_min)$beta
glm_fit <- glmnet(x, y, alpha = alpha_min)
plot(glm_fit)