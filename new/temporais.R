library(MASS)
library(bbmle)
source("plot-profmle.r")
temp <- read.csv("temporais.csv", sep=",", header=TRUE)
n <- function(t) temp[temp$time==t,3]
m <- function (t) log(n(t))
nLL = function(v1=1, v2=1, v3=1, c12=0.7, c13=0, c23=0, p1=-50, g1=2, p2=0, g2=0, p3=2, f=-10) {
	p1 <- exp(p1)/(1+exp(p1))
	p2 <- exp(p2)/(1+exp(p2))
	p3 <- exp(p3)/(1+exp(p3))
	g1 <- exp(g1)/(1+exp(g1))
	g2 <- exp(g2)/(1+exp(g2))
	f <- exp(f)
	v1 <- exp(v1)/(1+exp(v1)); v2 <- exp(v2)/(1+exp(v2)); v3 <- exp(v3)/(1+exp(v3))
	c12 <- -1 + 2 * exp(c12) / (1+exp(c12))
	c23 <- -1 + 2 * exp(c23) / (1+exp(c23))
	c13 <- -1 + 2 * exp(c13) / (1+exp(c13))
	hm <- function (t) log(c(p1*n(t)[1]+f*n(t)[3], g1*n(t)[1]+p2*n(t)[2], g2*n(t)[2] + p3*n(t)[3]))
	Sigma = matrix(c(v1, c12, c13, c12, v2, c23, c13, c23, v3), nrow=3)
	if (sum(is.nan((Sigma))) > 0) return (Inf)
	total = 0
	for (i in 1:5) total = total + (m(i) - hm(i)) %*% ginv(Sigma) %*% (m(i) - hm(i))
	5*log(det(Sigma)) + total
}

model <- mle2(nLL)
