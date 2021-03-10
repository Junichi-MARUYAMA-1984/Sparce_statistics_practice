# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

# Chapter2 e.g.11
# ロジスティック曲線の描画

f <- function(x) {
  return(exp(beta_0 + beta * x) / (1 + exp(beta_0 + beta * x)))
}
beta_0 <- 0
beta_seq <- c(0, 0.2, 0.5, 1, 2, 10)
m <- length(beta_seq)
beta <- beta_seq[1]
plot(f, xlim = c(-10, 10), ylim = c(0, 1),
     xlab = "x", ylab = "y",
     col = 1,
     main = "ロジスティック曲線")
for (i in 2:m) {
  beta <- beta_seq[i]
  par(new = TRUE)
  plot(f, xlim = c(-10, 10), ylim = c(0, 1),
       xlab = "", ylab = "",
       axes = FALSE, 
       col = i)
}
legend("topleft", legend = beta_seq,
       col = 1:m,
       lwd = 2,
       cex = 0.8)
par(new = FALSE)