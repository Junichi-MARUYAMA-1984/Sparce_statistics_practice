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
linear_lasso <- function(X, y, lambda = 0, beta = rep(0, ncol(X))) {
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
    print(eps)
    beta_old <- beta
  }
  
  beta <- beta / res$X_sd
  beta_0 <- res$y_bar - sum(res$X_bar * beta)
  
  return(list(beta = beta, beta_0 = beta_0))
}

# Function for Warm Start
warm_start <- function(X, y, lambda_max = 100) {
  dec <- round(lambda_max / 50)
  lambda_seq <- seq(lambda_max, 1, -dec)
  r <- length(lambda_seq)
  p <- ncol(X)
  
  coef_seq <- matrix(nrow = r, ncol = p)
  coef_seq[1, ] <- linear_lasso(X, y, lambda_seq[1])$beta
  
  for (k in 2:r) {
    coef_seq[k, ] <- linear_lasso(X, y, lambda_seq[k], coef_seq[(k - 1), ])$beta
  }
  
  return(coef_seq)
}

# Main routine

crime <- read.table("crime.txt")
X <- crime[, 3:7]
y <- crime[, 1]

coef_seq <- warm_start(X, y, 200)

p <- ncol(X)
lambda_max <- 200
dec <- round(lambda_max / 50)
lambda_seq <- seq(lambda_max, 1, -dec)

plot(log(lambda_seq), 
     coef_seq[, 1],
     ylim = c(min(coef_seq), max(coef_seq)), 
     xlab = "log(lambda)",
     ylab = "係数", 
     type = "n") 

for (j in 1:p) {
  lines(log(lambda_seq), coef_seq[, j], col = j)
}

legend("topright",
       legend = c("警察への年間資金",
                  "25歳以上で高校を卒業した人の割合",
                  "16~19歳で高校に通っていない人の割合", 
                  "18~24歳で大学生の割合", 
                  "25歳以上で4年制大学を卒業した人の割合"),
       col = 1:p,
       lwd = 2,
       cex = 0.8)
