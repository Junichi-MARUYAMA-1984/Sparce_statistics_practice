# Load libraries
library(survival)

# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

# Chapter2 e.g.20
# kidneyデータセットでの生存時間解析(1)
data(kidney)
y <- kidney$time
delta <- kidney$status
Surv(y, delta)