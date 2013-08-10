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

nLL = function(v1=0.05, v2=0.01, c12=0, p1=logit(0.3), g1=logit(0.3), p2=logit(0.4), f=log(1.45)) {
	p1 <- exp(p1)/(1+exp(p1))
	p2 <- exp(p2)/(1+exp(p2))
	g1 <- exp(g1)/(1+exp(g1))
	f <- exp(f)
	hm <- function (t) log(c(p1*n(t)[1]+f*n(t)[2], g1*n(t)[1]+p2*n(t)[2]))
	Sigma = matrix(c(v1, c12, c12, v2), nrow=2)
	total = 0
	for (i in 1:max(temp$time)) total = total + (m(i) - hm(i)) %*% ginv(Sigma) %*% (m(i) - hm(i))
	print (matrix(c(p1, f, g1, p2), ncol=2, byrow=T))
	return(as.numeric(max(temp$time)*log(det(Sigma)) + total))
}

#model <- mle2(nLL, optimizer="optimx", fixed=list(v1=0.05, v2=0.01, c12=0),
#			  lower=c(-5, -5, -5, -5), upper=c(5,5,5,5))
#model
#plot.profmle(profile(model))

vnLL <- Vectorize(nLL)
par(mfrow=c(2,2))
plot(logistic(seq(-1.5,0, by=0.1)), vnLL(p1=seq(-1.5,0,by=0.1)), type='l', main="p1")
plot(exp(seq(-1,1.2,by=0.1)),vnLL(f=seq(-1,1.2,by=0.1)), type='l', main="f")
plot(logistic(seq(-1.5,0, by=0.1)), vnLL(g1=seq(-1.5,0,by=0.1)), type='l', main="g1")
plot(logistic(seq(-1,0.4,by=0.1)),vnLL(p2=seq(-1,0.4,by=0.1)), type='l', main="p2")
