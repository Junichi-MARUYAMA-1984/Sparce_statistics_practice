# Load libraries
library(glmnet)

# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

# Chapter2 e.g.17
# Fisherのあやめデータを用いたMultiple logistic lasso解析
df <- iris
x <- as.matrix(df[, 1:4])
y <- as.vector(df[, 5])
n <- length(y)
u <- array(dim = n)
for (i in 1:n) {
  if (y[i] == "setosa") {
    u[i] <- 1
  } else if (y[i] == "versicolor") {
    u[i] <- 2
  } else {
    u[i] <- 3
  }
}
u <- as.numeric(u)
cv <- cv.glmnet(x, y, family = "multinomial")
cv2 <- cv.glmnet(x, y, family = "multinomial", type.measure = "class")
par(mfrow = c(1, 2))
plot(cv)
plot(cv2)
par(mfrow = c(1, 1))
lambda <- cv$lambda.min
result <- glmnet(x, y, lambda = lambda, family = "multinomial")
beta <- result$beta
beta_0 <- result$a0
v <- rep(0, n)
for (i in 1:n) {
  max_value <- -Inf
  for (j in 1:3) {
    value <- beta_0[j] + sum(beta[[j]] * x[i, ])
    if (value > max_value) {
      v[i] <- j
      max_value <- value
    }
  }
}
table(u, v)