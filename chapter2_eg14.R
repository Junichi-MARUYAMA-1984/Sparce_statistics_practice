# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

# Chapter2 e.g.14
# logistic lasso回帰による推定パラメタを用いての分類試行

# Function for determining formula type
soft_th <- function(lambda, x) {
  # xの絶対値がlambdaの絶対値よりも小さいときは、0を返す。
  # xの絶対値がlambdaの絶対値よりも大きいときは、
  # 二者の絶対値の差に「xの元々の符号」を付けて返す。
  return(sign(x) * pmax(abs(x) - lambda, 0))
}

# Function for centralization
centralize <- function(X, y, standardize = TRUE) {
  X <- as.matrix(X) # DataFrame型の引数Xを行列型に変換。
  n <- nrow(X) # 行列Xの行数
  p <- ncol(X) # 行列Xの列数
  X_bar <- array(dim = p) # 行列Xの各列成分の平均値を格納する配列
  X_sd <- array(dim = p) # 行列Xの各列成分の標準偏差を格納する配列
  
  for (j in 1:p) {
    X_bar[j] <- mean(X[, j]) # 行列X第j列成分の平均値を計算。
    X[, j] <- (X[, j] - X_bar[j]) # 行列X第j列の各成分に対して平均値を引き算。
    X_sd[j] <- sqrt(var(X[, j])) # 行列X第j列成分の標準偏差を計算。
    
    if (standardize == TRUE) {
      X[, j] <- X[, j] / X_sd[j] # 正規化
    }
  }
  
  if (class(y) == "matrix") { # yが行列の時は、各列に対して各列成分平均値を引き算。
    K <- ncol(y)
    y_bar <- array(dim = K)
    
    for (k in 1:K) {
      y_bar[k] <- mean(y[, k])
      y[, k] <- y[, k] - y_bar[k]
    }
  } else { # yがベクトルの時は、ベクトル成分の平均値を各成分から引き算。
    y_bar <- mean(y)
    y <- y - y_bar
  }
  
  return(list(X = X, y = y, X_bar = X_bar, X_sd = X_sd, y_bar = y_bar))
}

# Function for linear Lasso
linear_lasso <- function(X, y, lambda = 0, beta = rep(0, ncol(X))) {
  n <- nrow(X) # データフレームXの行数
  p <- ncol(X) # データフレームXの列数
  
  res <- centralize(X, y) # 引数X, yを中心化。
  X <- res$X # 中心化されたX。
  y <- res$y # 中心化されたy。
  
  eps <- 1 # betaの収束判定用変数を1で初期化。
  beta_old <- beta # beta_oldを引数betaで初期化。
  
  # Lassoのコアルーチン
  while (eps > 0.001) {
    for (j in 1:p) {
      r <- y - as.matrix(X[, -j]) %*% beta[-j]
      beta[j] <- soft_th(lambda, sum(r * X[, j]) / n) / (sum(X[, j] * X[, j]) / n)
    }
    eps <- max(abs(beta - beta_old))
    # print(eps)
    beta_old <- beta
  }
  
  beta <- beta / res$X_sd # 各betaの成分を正規化する前の値に戻す
  beta_0 <- res$y_bar - sum(res$X_bar * beta) # Lassoの切片beta_0
  
  return(list(beta = beta, beta_0 = beta_0))
}

# Weighted linear lasso
W_linear_lasso <- function(X, y, W, lambda = 0) {
  n <- nrow(X)
  p <- ncol(X)
  X_bar <- array(dim = p)
  for (k in 1:p) {
    X_bar[k] <- sum(W %*% X[, k]) / sum(W)
    X[, k] <- X[, k] - X_bar[k]
  }
  y_bar <- sum(W %*% y) / sum(W)
  y <- y - y_bar
  L <- chol(W) # W = t(L) %*% L
  u <- as.vector(L %*% y)
  V <- L %*% X
  beta <- linear_lasso(V, u, lambda)$beta
  beta_0 <- y_bar - sum(X_bar * beta)
  return(c(beta_0, beta))
}

# logistic lasso function
logistic_lasso <- function(X, y, lambda) {
  p <- ncol(X)
  beta <- Inf
  gamma <- rnorm(p)
  while (sum((beta - gamma) ^ 2) > 0.01) {
    beta <- gamma
    s <- as.vector(X %*% beta)
    v <- as.vector(exp(-s * y))
    u <- y * v / (1 + v)
    w <- v / (1 + v) ^ 2
    z <- s + u / w
    W <- diag(w)
    gamma <- W_linear_lasso(X[, 2:p], z, W, lambda = lambda)
    # print(gamma)
  }
  return(gamma)
}

# Main routine
# データ生成
N <- 100
p <- 2
X <- matrix(rnorm(N * p), ncol = p)
X <- cbind(rep(1, N), X)
beta <- 10 * rnorm(p + 1)
y <- array(N)
s <- as.vector(X %*% beta)
prob <- 1 / (1 + exp(s))
for (i in 1:N) {
  if (runif(1) > prob[i]) {
    y[i] <- 1
  } else {
    y[i] <- -1
  }
}

# パラメタ推定
beta_est <- logistic_lasso(X, y, 0.1)

# 分類処理
for (i in 1:N) {
  if (runif(1) > prob[i]) {
    y[i] <- 1
  } else {
    y[i] <- -1
  }
}
z <- sign(X %*% beta_est) # 指数部が正ならy = 1, 負ならy = -1と推定する
table(y, z)