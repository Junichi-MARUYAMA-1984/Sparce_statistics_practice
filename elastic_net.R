# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

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

# Function for elastic net
elastic_net <- function(X, y, lambda = 0, beta = rep(0, ncol(X)), alpha = 1) {
  n <- nrow(X) # データフレームXの行数
  p <- ncol(X) # データフレームXの列数
  
  res <- centralize(X, y) # 引数X, yを中心化。
  X <- res$X # 中心化されたX。
  y <- res$y # 中心化されたy。
  
  eps <- 1 # betaの収束判定用変数を1で初期化。
  beta_old <- beta # beta_oldを引数betaで初期化。
  
  # Elastic netのコアルーチン
  while (eps > 0.001) {
    for (j in 1:p) {
      r <- y - as.matrix(X[, -j]) %*% beta[-j]
      beta[j] <- soft_th(lambda * alpha, 
                         sum(r * X[, j]) / n) / (sum(X[, j] * X[, j]) / n + lambda * (1 - alpha))
    }
    eps <- max(abs(beta - beta_old))
    print(eps)
    beta_old <- beta
  }
  
  beta <- beta / res$X_sd # 各betaの成分を正規化する前の値に戻す
  beta_0 <- res$y_bar - sum(res$X_bar * beta) # Lassoの切片beta_0
  
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
  coef_seq[i, ] <- elastic_net(X, y, lambda_seq[i], alpha = 0)$beta
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
