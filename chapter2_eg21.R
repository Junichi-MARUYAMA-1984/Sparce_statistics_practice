# Load libraries
library(survival)

# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

# Chapter2 e.g.20
# kidneyデータセットでの生存時間解析(2)
fit <- survfit(Surv(time, status) ~ disease, data = kidney)
plot(fit, xlab = "時間", ylab = "生存率",
     col = c("red", "green", "blue", "black"))
legend(300, 0.8,
       legend = c("その他", "GN", "AN", "PKD"),
       lty = 1,
       col = c("red", "green", "blue", "black"))