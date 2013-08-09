library(MASS)
library(bbmle)
library(optimx)
source("plot-profmle.r")
source("temporais-simula.R")
#temp <- read.csv("temporais.csv", sep=",", header=TRUE)
n <- function(t) temp[temp$time==t,3]
m <- function (t) log(n(t))

logistic <- function (p) exp(p)/(1+exp(p))
logit <- function (p) log(p) - log(1-p)

nLL = function(v1=1, v2=1, v3=1, c12=0.7, c13=0, c23=0, p1=-50, g1=2, p2=0, g2=0, p3=2, f=-10) {
	p1 <- exp(p1)/(1+exp(p1))
	p2 <- exp(p2)/(1+exp(p2))
	p3 <- exp(p3)/(1+exp(p3))
	g1 <- exp(g1)/(1+exp(g1))
	g2 <- exp(g2)/(1+exp(g2))
	f <- exp(f)
	hm <- function (t) log(c(p1*n(t)[1]+f*n(t)[3], g1*n(t)[1]+p2*n(t)[2], g2*n(t)[2] + p3*n(t)[3]))
	Sigma = matrix(c(v1, c12, c13, c12, v2, c23, c13, c23, v3), nrow=3)
	total = 0
	for (i in 1:max(temp$time)) total = total + (m(i) - hm(i)) %*% ginv(Sigma) %*% (m(i) - hm(i))
	print (matrix(c(p1, 0, f, g1, p2, 0, 0, g2, p3), ncol=3, byrow=T))
	return(as.numeric(5*log(det(Sigma)) + total))
}

model <- mle2(nLL, optimizer="optimx",
			  start=list(p1=logit(0.3), g1=logit(0.3), f=log(1.5), g2=logit(0.4), p2=logit(0.4), p3=logit(0.6)),
			  	   fixed=list(v1=0.1,v2=0.1,v3=0.1,
							  c12=0,c13=0,c23=0
							  ))
model
plot.profmle(profile(model))
