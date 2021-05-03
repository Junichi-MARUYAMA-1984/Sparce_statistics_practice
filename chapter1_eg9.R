# Load libraries
library(glmnet)

# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

df <- read.table("crime.txt")
X <- as.matrix(df[, 3:7])
y <- as.vector(df[, 1])
cv_fit <- cv.glmnet(X, y)
plot(cv_fit)
lambda_min <- cv_fit$lambda.min
lambda_min