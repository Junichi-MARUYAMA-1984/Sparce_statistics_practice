# Load libraries
library(MASS)

# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

# Function for the squared coefficient of determination (R^2)
R2 <- function(x, y) {
  y_hat <- lm(y ~ x)$fitted.values
  y_bar <- mean(y)
  
  RSS <- sum((y - y_hat)^2)
  TSS <- sum((y - y_bar)^2)
  
  return(1 - RSS / TSS)
}

# Function for Variance Inflation Factor
vif <- function(x) {
  p <- ncol(x)
  values <- array(dim = p)
  
  for (j in 1:p) {
    values[j] <- 1 / (1 - R2(x[, -j], x[, j]))
  }
  
  return(values)
}

# Main routine
x <- as.matrix(Boston)
vif(x)
