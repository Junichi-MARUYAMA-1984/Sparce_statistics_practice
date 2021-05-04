# Load libraries
library(glmnet)
library(MASS)

# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

# Chapter2 e.g.19
# birthwtデータを用いたpoisson lasso解析
data(birthwt)
df <- birthwt[, -1]
dy <- df[, 8]
dx <- data.matrix(df[, -8])
cvfit <- cv.glmnet(x = dx, y = dy, family = "poisson", standardize = TRUE)
coef(cvfit, s = "lambda.min")