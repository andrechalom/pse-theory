library(MASS)
#library(bbmle)
#library(optimx)
#source("plot-profmle.r")
#source("temporais-simula.R")
temp <- read.csv("temporais.csv", sep=",", header=TRUE)
n <- function(t) temp[temp$time==t,3]
m <- function (t) log(n(t))

logistic <- function (p) exp(p)/(1+exp(p))
logit <- function (p) log(p) - log(1-p)

nLL = function(v1=0.03, v2=0.01, v3=0.01, c12=-0.013, c13=-0.014, c23=0.009, p1=2.67, g1=-2.67, p2=3.41, g2=-1.42, p3=3.79, f=-7.93) {
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
	return(as.numeric(max(temp$time)*log(det(Sigma)) + total))
}
print(nLL())

#model <- mle2(nLL, optimizer="optimx",
#			  start=list(p1=logit(0.3), g1=logit(0.3), f=log(1.5), g2=logit(0.4), p2=logit(0.4), p3=logit(0.6)),
#			  	   fixed=list(v1=0.1,v2=0.1,v3=0.1,
#							  c12=0,c13=0,c23=0
#							  ))
#model
#plot.profmle(profile(model))
inc=0.01
vnLL <- Vectorize(nLL)
par(mfrow=c(4,3))
sq = seq(-5,5,by=inc)
p1 = vnLL(p1=sq)
print(sq[which(p1==min(p1))])
plot(logistic(sq), p1, type='l', main="p1")
sq = seq(3,5,by=inc)
p2 = vnLL(p2=sq)
print(sq[which(p2==min(p2))])
plot(logistic(sq), p2, type='l', main="p2")
sq = seq(-5,5,by=inc)
p3 = vnLL(p3=sq)
print(sq[which(p3==min(p3))])
plot(logistic(sq), p3, type='l', main="p3")
sq = seq(-5,5,by=inc)
g1 = vnLL(g1=sq)
print(sq[which(g1==min(g1))])
plot(logistic(sq), g1, type='l', main="g1")
sq = seq(-5,5,by=inc)
g2 = vnLL(g2=sq)
print(sq[which(g2==min(g2))])
plot(logistic(sq), g2, type='l', main="g2")
sq = seq(-8,-4, by=inc)
f = vnLL(f=sq)
print(sq[which(f==min(f))])
plot(exp(sq), f, type='l', main="f")
sq = seq(inc,2, by=inc)
v1 = vnLL(v1=sq)
print(sq[which(v1==min(v1[!is.nan(v1)]))])
plot(sq, v1, type='l', main="v1")
v2 = vnLL(v2=sq)
print(sq[which(v2==min(v2))])
plot(sq, v2, type='l', main="v2")
v3 = vnLL(v3=sq)
print(sq[which(v3==min(v3))])
plot(sq, v3, type='l', main="v3")
sq = seq(0,60, by=inc)/1000
c12 = vnLL(c12=sq)
print(sq[which(c12==min(c12[!is.nan(c12)]))])
plot(sq, c12, type='l', main="c12")
c13 = vnLL(c13=sq)
print(sq[which(c13==min(c13[!is.nan(c13)]))])
plot(sq, c13, type='l', main="c13")
c23 = vnLL(c23=sq)
print(sq[which(c23==min(c23[!is.nan(c23)]))])
plot(sq, c23, type='l', main="c23")

# Simulando a populacao
p1=logistic(2.67); g1=logistic(-2.67); p2=logistic(3.41); g2=logistic(-3.42); p3=logistic(3.79); f=exp(-7.93);
m <- matrix( c(p1, 0, f, g1, p2, 0, 0, g2, p3), nrow=3, byrow=T)
n0 = temp[temp$time==0,3]
simul <- data.frame(time=c(0,0,0), stage=c(1,2,3), density=n0)
for (i in 1:5) {
	n0 = m%*%n0
	simul <- rbind(simul, c(i, 1, n0[1]))
	simul <- rbind(simul, c(i, 2, n0[2]))
	simul <- rbind(simul, c(i, 3, n0[3]))
}

plot(temp[temp$stage==1,3], lty=1, lwd=2, type='l', ylim=c(0, 170))
points(temp[temp$stage==2,3], lty=1, lwd=2, type='l', col='red')
points(temp[temp$stage==3,3], lty=1, lwd=2, type='l', col='green')
points(simul[simul$stage==1,3], lty=3, type='l')
points(simul[simul$stage==2,3], lty=3, type='l', col='red')
points(simul[simul$stage==3,3], lty=3, type='l', col='green')
