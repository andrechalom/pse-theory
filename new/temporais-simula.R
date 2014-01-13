##### Matriz com 3x3: reverter para versao 1066
library(MASS)
library(expm)
logistic <- function (p) exp(p)/(1+exp(p))
logit <- function (p) log(p) - log(1-p)
### Matriz "correta"
m <- matrix(c(0.3, 1.45, 0.3, 0.4), ncol=2, byrow=T)
correct <- c(logit(0.3), log(1.45), logit(0.3), logit(0.4))
#### | P1 | F  |
#### | G1 | P2 |
# Simulando uma serie temporal
n0 <- c(90, 43)
temp <- data.frame(time=c(0,0), stage=c(1,2), density=n0)
for (i in 1:200) {
	D <- matrix(c(rnorm(1, 0, 0.05), 0, 0, rnorm(1,0,0.1)), nrow=2, byrow=T)
	n0 <- expm(D)%*%m%*%n0
	temp <- rbind(temp, c(i, 1, n0[1]))
	temp <- rbind(temp, c(i, 2, n0[2]))
}
plot(temp[temp$stage==1,3])
x <- data.frame(temp[temp$stage==1,3], temp[temp$stage==2,3])
matplot(x)

Simula <- function(p1, f, g1, p2) {
		Np <- c(90, 43) #Pop inicial!
		myp <- Np
		for(i in 1:200) {
				L <- matrix(c(p1, f, g1, p2), nrow=2, byrow=TRUE)
				myp <- L %*% myp
				Np <- c(Np, myp)
		}
	return (Np);
}

### OTIMIZACAO NAO FUNCIONA!!!
constrain <- function(x) c(logistic(x[1]), exp(x[2]), logistic(x[3]), logistic(x[4]))
ssqres <- function (x) sum((Simula(x[1], x[2], x[3], x[4]) - temp[,3])^2)
cssqres <- function(x) ssqres(constrain(x))
a <- optim(par=c(-1, 0.2, -1, 0), fn=cssqres, method="SANN", control=list(maxit=200000))
a <- optim(a$par, fn=cssqres, method="SANN", control=list(maxit=2000000))
#a <- c()
#a$par <- c(-3, -14.8, -3.8, -12)
test <- Simula(exp(a$par[1]), logistic(a$par[2]), logistic(a$par[3]), logistic(a$par[4]))
test <- matrix(test, ncol=2, byrow=T)
matplot(test)
#minsq <- sum((test-temp[,3])^2)
#sum((test-temp[,3])^2)-minsq
