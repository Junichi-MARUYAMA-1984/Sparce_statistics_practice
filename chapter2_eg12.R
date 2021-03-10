# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

# Chapter2 e.g.12
# ロジスティック回帰パラメータの最尤推定

# データ生成
N <- 1000
p <- 2
X <- matrix(rnorm(N * p), ncol = p)
X <- cbind(rep(1, N), X)
beta <- rnorm(p + 1)
y <- array(N)
s <- as.vector(X %*% beta)
prob <- 1 / (1 + exp(s)) # 式(2.4)で定義されるP(Y = -1|x)
# 0 < rand < 1なる一様乱数randを発生させた時、
# rand > probとなる確率は、P(Y = 1|x)に等しい。
for (i in 1:N) {
  if (runif(1) > prob[i]) {
    y[i] <- 1
  } else {
    y[i] <- -1
  }
}
print(beta)

# 最尤推定値の計算
beta_est <- Inf
gamma <- rnorm(p + 1)
while (sum((beta_est - gamma) ^ 2) > 0.001) {
  beta_est <- gamma
  s <- as.vector(X %*% beta_est)
  v <- exp(-s * y)
  u <- y * v / (1 + v)
  w <- v / (1 + v) ^ 2
  z <- s + u / w
  W <- diag(w)
  gamma <- as.vector(solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% z)
  print(gamma)
}
print(beta)
