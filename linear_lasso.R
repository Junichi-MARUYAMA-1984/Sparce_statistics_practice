# Function for determining formula type
soft_th <- function(lambda, x) {
  return(sign(x) * pmax(abs(x) - lambda, 0))
}

# Function for centralization
centralize <- function(X, y, standardize = TRUE) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  X_bar <- array(dim = p)
  X_sd <- array(dim = p)
  
  for (j in 1:p) {
    X_bar[j] <- mean(X[, j])
    X[, j] <- (X[, j] - X_bar[j])
    X_sd[j] <- sqrt(var(X[, j]))
    
    if (standardize == TRUE) {
      X[, j] <- X[, j] / X_sd[j]
    }
  }
  
  if (class(y) == "matrix") { # in case y is matrix
    K <- ncol(y)
    y_bar <- array(dim = K)
    
    for (k in 1:K) {
      y_bar[k] <- mean(y[, k])
      y[, k] <- y[, k] - y_bar[k]
    }
  } else { # in case y is vector
    y_bar <- mean(y)
    y <- y - y_bar
  }
  
  return(list(X = x, y = y, X_bar = X_bar, X_sd = X_sd, y_bar = y_bar))
}

# Function for linear Lasso
linear_lasso <- function(X, y, lambda, beta) {
  lambda <- 0
  beta <- rep(0, ncol(X))
  n <- nrow(X)
  p <- ncol(X)
  
  res <- centralize(X, y)
  X <- res$X
  y <- res$y
  
  eps <- 1
  beta_old <- beta
  
  while (eps > 0.001) {
    for (j in 1:p) {
      r <- y - as.matrix(X[, -j]) %*% beta[-j]
      beta[j] <- soft_th(lambda, sum(r * X[, j]) / n) / (sum(X[, j] * X[, j]) / n)
    }
    eps <- max(abs(beta - beta_old))
    beta_old <- beta
  }
  
  beta <- beta / res$X_sd
  beta_0 <- res$y_bar - sum 
}