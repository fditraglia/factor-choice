setwd("~/factor-choice/")
library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
sourceCpp("sampler.cpp")
source("calibrated_simulation.R")

value <- read.csv("data_value.csv")
equal <- read.csv("data_equal.csv")
FF3 <- c("Mkt.RF", "SMB", "HML")
Size10 <- c("Lo10", paste0("Dec", 2:9), "Hi10")
Industry10 <- c("NoDur", "Durbl", "Manuf", "Enrgy", 
                "HiTec", "Telcm", "Shops", "Hlth", 
                "Utils", "Other")

set.seed(8372)

G0 <- 10 * diag(40)
g0 <- rep(0, 40)
R0 <- diag(10)
r0 <- 15

valueSize <- SURidentical(FF3, Size10, value)
Gamma <- valueSize$Gamma
Sigma <- valueSize$Sigma
errors <- rmvnorm(nrow(value), sigma = Sigma)
X <- cbind(alpha = rep(1, nrow(value)), as.matrix(value[,FF3]))
Y <- X %*% Gamma + errors

system.time(gibbsValueSize <- samplerTest(X, Y, G0, g0, R0, r0, 4000, 1000))
g_star <- matrix(rowMeans(gibbsValueSize$g_draws), 4, 10)
Omega_inv_star <- rowMeans(gibbsValueSize$Omega_inv_draws)

g_star - Gamma
solve(devech(Omega_inv_star, 10)) - Sigma

