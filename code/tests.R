library("Rcpp")
library("RcppArmadillo")
sourceCpp("distributions.cpp")

library("MCMCpack")
library("mvtnorm")

M <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
M_inv <- matrix(c(4/3, -2/3, -2/3, 4/3), 2, 2)
m <- c(-3, 3)

set.seed(342)
foo <- replicate(100, rwish(5, M))
set.seed(342)
bar <- replicate(100, draw_wishart(5, M))
all.equal(foo, bar)

set.seed(1234)
foo <- t(drop(replicate(10000, draw_normal(m, M_inv))))
cov(foo) - M
colMeans(foo) - m
qqnorm(rowSums(foo))

x <- cbind(c(1, 1), c(-1, 1), c(-1, -1))
foo <- density_normal(x, mu = m, Sigma_inv = M_inv)
bar <- dmvnorm(t(x), mean = m, sigma = M)
all.equal(drop(foo), bar)
boo <- density_normal(x, mu = m, Sigma_inv = M_inv, logret = TRUE)
baz <- dmvnorm(t(x), mean = m, sigma = M, log = TRUE)
all.equal(drop(boo), baz)

foo <- density_wishart(3 * M_inv, 5, M_inv)
bar <- dwish(3 * M_inv, 5, M_inv) 
all.equal(foo, bar)
baz <- density_wishart(3 * M_inv, 5, M_inv, logret = TRUE)
all.equal(baz, log(foo))
