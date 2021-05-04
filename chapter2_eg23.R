# Load libraries
library(glmnet)
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

# Clear Working Memory
rm(list = ls())

# Setting Japanese font (for MacOSX)
# par(family= "HiraKakuProN-W3")

# Chapter2 e.g.23
# 悪性リンパ腫データを用いたcox lasso解析

load("LymphomaData.rda")
attach("LymphomaData.rda")
names(patient.data)
x <- t(patient.data$x)
y <- patient.data$time
delta <- patient.data$status
Surv(y, delta)
cv_fit <- cv.glmnet(x, Surv(y, delta), family = "cox")
fit2 <- glmnet(x, Surv(y, delta), lambda = cv_fit$lambda.min, family = "cox")
z <- sign(drop(x %*% fit2$beta))
fit3 <- survfit(Surv(y, delta) ~ z)
autoplot(fit3)
mean(y[z == 1])
mean(y[z == -1])