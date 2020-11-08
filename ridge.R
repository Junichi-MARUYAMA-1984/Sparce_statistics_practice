# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

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
  
  return(list(X = X, y = y, X_bar = X_bar, X_sd = X_sd, y_bar = y_bar))
}

# Function for least squares method
linear <- function(X, y) {
  X <- as.matrix(X)
  p <- ncol(X) # 行列Xの列数
  n <- nrow(X) # 行列Xの行数
  
  res <- centralize(X, y) # 引数X, yを中心化。
  X <- res$X # 中心化されたX
  y <- res$y # 中心化されたy
  
  # LSMのコア処理
  beta <- as.vector(solve(t(X) %*% X) %*% t(X) %*% y)
  
  beta <- beta / res$X_sd # 各betaの成分を正規化する前の値に戻す。
  beta_0 <- res$y_bar - sum(res$X_bar * beta) # LSMの切片beta_0
  
  return(list(beta = beta, beta_0 = beta_0))
}

# Function for Ridge
ridge <- function(X, y, lambda = 0) {
  X <- as.matrix(X)
  p <- ncol(X) # 行列Xの列数
  n <- nrow(X) # 行列Xの行数
  
  res <- centralize(X, y) # 引数X, yを中心化。
  X <- res$X # 中心化されたX
  y <- res$y # 中心化されたy
  
  # Ridgeのコア処理
  beta <- drop(solve(t(X) %*% X + n * lambda * diag(p)) %*% t(X) %*% y)
  
  beta <- beta / res$X_sd # 各betaの成分を正規化する前の値に戻す。
  beta_0 <- res$y_bar - sum(res$X_bar * beta) # Ridgeの切片beta_0
  
  return(list(beta = beta, beta_0 = beta_0))
}

# Main routine
df <- read.table("crime.txt")
X <- df[, 3:7]
y <- df[, 1]

p <- ncol(X)
lambda_seq <- seq(0, 200, 0.1)

plot(lambda_seq, 
     xlim = c(0, 200),
     ylim = c(-10, 20), 
     xlab = "lambda",
     ylab = "beta", 
     main = "各lambdaにおける各項目beta値\n（年間犯罪率への寄与率を表す）",
     type = "n",
     col = "red")

r <- length(lambda_seq)
coef_seq <- array(dim = c(r, p)) # 各lambdaに対するbeta値を格納する配列

for (i in 1:r) {
  coef_seq[i, ] <- ridge(X, y, lambda_seq[i])$beta
}

for (j in 1:p) {
  par(new = TRUE)
  lines(lambda_seq, coef_seq[, j], col = j)
}

legend("topright",
       legend = c("警察への年間資金",
                  "25歳以上で高校を卒業した人の割合",
                  "16~19歳で高校に通っていない人の割合", 
                  "18~24歳で大学生の割合", 
                  "25歳以上で4年制大学を卒業した人の割合"),
       col = 1:p, # 凡例中の線の色（プロットされている線の色と揃える）
       lwd = 2, # 凡例中の線の太さ
       cex = 0.8) # 凡例フォントの大きさ
