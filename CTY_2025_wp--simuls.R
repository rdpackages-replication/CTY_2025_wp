rm(list=ls(all=TRUE))

library(MASS)
library(tidyr)
library("lmtest")
library("sandwich")
library(ggplot2)
library(rdrobust)
library(latex2exp)
library(sf)
library(dplyr)  # For sampling
library(progressr)
library(doParallel)
library(foreach)
library(expm)

set.seed(3)

################################## Load Data ###################################

x <- read.csv("CIT_2023_CUP_multiscore-nongeo.csv")
x$X <- NULL
colnames(x) <- c("x.1", "x.2","y","d")
na.ok <- complete.cases(x$x.1) & complete.cases(x$x.2)
x <- x[na.ok,]

neval <- 40
eval <- matrix(nrow = neval, ncol = 2)
for (i in 1: ceiling(neval * 0.5)){
  eval[i,] <- c(0, 50 - (i-1) * 50 / ceiling(neval * 0.5))
}
for (i in (ceiling(neval * 0.5)+1): neval){
  eval[i,] <- c((i - ceiling(neval * 0.5) - 1) *50 / (ceiling(neval * 0.5)),0)
}
eval <- data.frame(eval)
colnames(eval) <- c("x.1", "x.2")

######################## Polynomial Fits to the Data ###########################

formula <- y ~ x.1 + x.2 # run for DGP 1 
# formula <- y ~ x.1 + x.2 + I(x.1^2) + I(x.1 * x.2) + I(x.2^2) # run for DGP 2

# Fit models separately for d = 0 and d = 1
model_d0 <- lm(formula, data = x[x$d == 0, ])
model_d1 <- lm(formula, data = x[x$d == 1, ])

# Summary of the models
summary(model_d0)
summary(model_d1)

# Simulate data with independent Beta-distributed coordinates
m <- 1000

alpha_1 <- 3
beta_1 <- 4
alpha_2 <- 3
beta_2 <- 4

n <- 20 * 20 * 25 * 2
# n <- 40 * 40 * 25
# n <- 10 * 10 * 25
h_grid <- c(25)

num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)  # Register parallel backend

system.time({
  foreach(j = 1:m) %dopar% {
    library(MASS)
    library(tidyr)
    library("lmtest")
    library("sandwich")
    library(ggplot2)
    library(rdrobust)
    library(latex2exp)
    library(sf)
    library(dplyr)  # For sampling
    library(progressr)
    library(doParallel)
    library(foreach)
    library(expm)
    library(rd2d)

    eval.subset <- eval
    est.2d <- array(0,dim = c(neval,1))
    est.1d <- array(0,dim = c(neval,1))
    ci.upper <- array(0,dim = c(neval,1))
    ci.lower <- array(0,dim = c(neval,1))
    cb.upper <- array(0,dim = c(neval,1))
    cb.lower <- array(0,dim = c(neval,1))

    # generate data set for the ith monte carlo

    x.1 <- 100 * rbeta(n, alpha_1, beta_1) - 25
    x.2 <- 100 * rbeta(n, alpha_2, beta_2) - 25

    X <- cbind(x.1, x.2)
    X <- as.data.frame(X)
    t <- x.1 >= 0 & x.2 >= 0
    y <- ifelse(t == 0, predict(model_d0, newdata = X), predict(model_d1, newdata = X)) * 2
    y <- y + rnorm(n,mean = 0, sd = 0.3306) * (1 - t) + rnorm(n,mean = 0, sd = 0.4348) * t

    D <- proxy::dist(X, eval.subset, method = "euclidean")  # Use "euclidean" for Euclidean distances
    d_expanded <- matrix(rep(2 * t - 1, times = ncol(D)), nrow = nrow(D), ncol = ncol(D))
    D <- D * d_expanded

    result.rd2d <- rd2d(y, X, t, eval.subset, stdvars = FALSE, repp = 5000, vce = "hc3")
    out.rd2d <- cbind(result.rd2d$opt$h01, result.rd2d$tau.hat, result.rd2d$cb$CI.l, result.rd2d$cb$CI.r, 
                      result.rd2d$cb$CB.l, result.rd2d$cb$CB.r, result.rd2d$tau.hat.q)

    result.dist.kinkon <- rd2d.dist(y,D,kink = "on", repp = 5000, vce = "hc3")
    out.dist.kinkon <- cbind(result.dist.kinkon$opt$h0, result.dist.kinkon$tau.hat, result.dist.kinkon$cb$CI.l, result.dist.kinkon$cb$CI.r,
                             result.dist.kinkon$cb$CB.l, result.dist.kinkon$cb$CB.r, result.dist.kinkon$tau.hat.q)

    result.dist.kinkoff <- rd2d.dist(y,D,kink = "off", repp = 5000, vce = "hc3")
    out.dist.kinkoff <- cbind(result.dist.kinkoff$opt$h0, result.dist.kinkoff$tau.hat, result.dist.kinkoff$cb$CI.l, result.dist.kinkoff$cb$CI.r,
                              result.dist.kinkoff$cb$CB.l, result.dist.kinkoff$cb$CB.r, result.dist.kinkoff$tau.hat.q)

    write.csv(out.rd2d, file = sprintf("rd2d_monte/rd2d_sim%d_m1000_sd01_n20000_linear.csv", j), row.names = FALSE) # run for DGP 1
    write.csv(out.dist.kinkon, file = sprintf("rd2d_monte/dist_kinkon_sim%d_m1000_sd01_n20000_linear.csv", j), row.names = FALSE) # run for DGP 1
    write.csv(out.dist.kinkoff, file = sprintf("rd2d_monte/dist_kinkoff_sim%d_m1000_sd01_n20000_linear.csv", j), row.names = FALSE) # run for DGP 1
    
    # write.csv(out.rd2d, file = sprintf("rd2d_monte/rd2d_sim%d_m1000_sd01_n20000_hc3.csv", j), row.names = FALSE) # run for DGP 2
    # write.csv(out.dist.kinkon, file = sprintf("rd2d_monte/dist_kinkon_sim%d_m1000_sd01_n20000_hc3.csv", j), row.names = FALSE) # run for DGP 2
    # write.csv(out.dist.kinkoff, file = sprintf("rd2d_monte/dist_kinkoff_sim%d_m1000_sd01_n20000_hc3.csv", j), row.names = FALSE) # run for DGP 2
  }})

stopCluster(cl)
