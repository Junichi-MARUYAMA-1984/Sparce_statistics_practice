# Load libraries
library(glmnet)

# Clear Working Memory
rm(list = ls())

# Main routine
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

glm_fit_lasso <- glmnet(x, y)
glm_fit_ridge <- glmnet(x, y, alpha = 0)

plot(glm_fit_lasso, xvar = "lambda")
legend("topleft",
       legend = c("X1", "X2", "X3", "X4", "X5", "X6"),
       col = 1:6,
       lwd = 2,
       cex = 0.8)

plot(glm_fit_ridge, xvar = "lambda")
legend("topleft",
       legend = c("X1", "X2", "X3", "X4", "X5", "X6"),
       col = 1:6,
       lwd = 2,
       cex = 0.8)