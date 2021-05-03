# Load libraries
library(glmnet)

# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

# Chapter2 e.g.15
# 乳がん遺伝子発現データを用いたlogistic lasso解析
df <- read.csv("breastcancer.csv")
x <- as.matrix(df[, 1:1000])
y <- as.vector(df[, 1001])
cv <- cv.glmnet(x, y, family = "binomial")
cv2 <- cv.glmnet(x, y, family = "binomial", type.measure = "class")
par(mfrow = c(1, 2))
plot(cv)
plot(cv2)
par(mfrow = c(1, 1))
glm <- glmnet(x, y, lambda = 0.03, family = "binomial")
beta <- drop(glm$beta)
beta[beta != 0]