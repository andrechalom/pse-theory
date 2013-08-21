library(MASS)
library(bbmle)
library(optimx)
library(mvtnorm)
#source("plot-profmle.r")
#source("temporais-simula.R")
temp <- read.csv("temporais.csv", sep=",", header=TRUE)
n <- function(t) temp[temp$time==t,3]
m <- function (t) log(n(t))

logistic <- function (p) exp(p)/(1+exp(p))
logit <- function (p) log(p) - log(1-p)

# nao foram: p1 p2 g1 f v1 v2 c12 c13
nLL = function(p1=4, p2=4, p3=5.27, g1=-1, g2=-2.48, f=-6, v1=-1, v2=-2, v3=-5.28, c12=0, c13=-0.014, c23=-6.18) {
	p1 <- exp(p1)/(1+exp(p1))
	p2 <- exp(p2)/(1+exp(p2))
	p3 <- exp(p3)/(1+exp(p3))
	g1 <- exp(g1)/(1+exp(g1))
	g2 <- exp(g2)/(1+exp(g2))
	f <- exp(f)
	v1 <- exp(v1); v2<- exp(v2); v3 <- exp(v3)
	c12 <- -1+2*logistic(c12)
	c23 <- -1+2*logistic(c23)
	c13 <- -1+2*logistic(c13)
	hm <- function (t) log(c(p1*(1-g1)*n(t)[1]+f*n(t)[3], p1*g1*n(t)[1]+p2*(1-g2)*n(t)[2], p2*g2*n(t)[2] + p3*n(t)[3]))
	Sigma = matrix(c(v1*v1, c12*v1*v2, c13*v1*v3, c12*v1*v2, v2*v2, c23*v2*v3, c13*v1*v3, c23*v2*v3, v3*v3), nrow=3)
	if(is.nan(log(det(Sigma))) | is.na(log(det(Sigma)))) return (0);
	total = 0
	for (i in 1:max(temp$time)) total = total + (m(i) - hm(i)) %*% ginv(Sigma) %*% (m(i) - hm(i))
	print(as.numeric(max(temp$time)*log(det(Sigma)) + total))
	print(dmvnorm(x=m(1), mean=hm(1), sigma=Sigma))
}

# model <- mle2(nLL, optimizer="optimx")
#               start=list(p1=logit(0.3), g1=logit(0.3), f=log(1.5), g2=logit(0.4), p2=logit(0.4), p3=logit(0.6),
#                          v1=0),
#                 fixed=list(v1=0.1,v2=0.1,v3=0.1,
#                      c12=0,c13=0,c23=0
#                               ))
# model
#plot.profmle(profile(model))
inc=0.01
vnLL <- Vectorize(nLL)
par(mfrow=c(4,3))
sq = seq(-7,7,by=inc)
p1 = vnLL(p1=sq)
plot(logistic(sq), p1, type='l', main="p1")
lp1 = min(p1)
p1 =sq[which(p1==min(p1))]
p2 = vnLL(p2=sq)
plot(logistic(sq), p2, type='l', main="p2")
lp2=min(p2)
p2=sq[which(p2==min(p2))]
p3 = vnLL(p3=sq)
plot(logistic(sq), p3, type='l', main="p3")
lp3=min(p3)
p3=sq[which(p3==min(p3))]
g1 = vnLL(g1=sq)
plot(logistic(sq), g1, type='l', main="g1")
lg1=min(g1)
g1=sq[which(g1==min(g1))]
g2 = vnLL(g2=sq)
plot(logistic(sq), g2, type='l', main="g2")
lg2=min(g2)
g2=sq[which(g2==min(g2))]
sq = seq(-8,-4, by=inc)
f = vnLL(f=sq)
plot(exp(sq), f, type='l', main="f")
lf = min(f)
f=sq[which(f==min(f))]
sq = seq(-7,7, by=inc)
v1 = vnLL(v1=sq)
plot(sq, v1, type='l', main="v1")
lv1=min(v1)
v1=sq[which(v1==min(v1[!is.nan(v1)]))]
v2 = vnLL(v2=sq)
plot(sq, v2, type='l', main="v2")
lv2=min(v2)
v2=sq[which(v2==min(v2))]
v3 = vnLL(v3=sq)
plot(sq, v3, type='l', main="v3")
lv3=min(v3)
v3=sq[which(v3==min(v3))]
c12 = vnLL(c12=sq)
plot(sq, c12, type='l', main="c12")
lc12=min(c12)
c12=sq[which(c12==min(c12[!is.nan(c12)]))]
c13 = vnLL(c13=sq)
plot(sq, c13, type='l', main="c13")
lc13=min(c13)
c13=sq[which(c13==min(c13[!is.nan(c13)]))]
c23 = vnLL(c23=sq)
plot(sq, c23, type='l', main="c23")
lc23=min(c23)
c23=sq[which(c23==min(c23[!is.nan(c23)]))]
print(data.frame(par=c("p1", "p2", "p3", "g1", "g2", "f", "v1", "v2", "v3", "c12", "c13", "c23"),
				 logl=round(nLL()-c(lp1, lp2, lp3, lg1, lg2, lf, lv1, lv2, lv3, lc12, lc13, lc23)),
				 value=c(p1, p2, p3, g1, g2, f, v1, v2, v3, c12, c13, c23)))


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
