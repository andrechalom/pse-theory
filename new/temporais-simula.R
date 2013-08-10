##### Matriz com 3x3: reverter para versao 1066
library(MASS)
library(bbmle)
library(optimx)
source("plot-profmle.r")
library(expm)
### Matriz "correta"
m <- matrix(c(0.3, 1.45, 0.3, 0.4), ncol=2, byrow=T)
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
