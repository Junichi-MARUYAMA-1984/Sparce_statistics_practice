# Load libraries
library(glmnet)
library(MASS)

# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

# Main routine
df <- Boston
x <- as.matrix(df[, 1:13])
y <- df[, 14]

fit <- glmnet(x, y)
plot(fit, xvar = "lambda", main = "BOSTON")