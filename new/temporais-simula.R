library(MASS)
library(bbmle)
library(optimx)
source("plot-profmle.r")
### Matriz "correta"
m <- matrix(c(0.3, 0, 1.5, 0.3, 0.4, 0, 0, 0.4, 0.6), ncol=3, byrow=T)
# Simulando uma serie temporal
n0 <- c(22, 53, 81)
temp <- data.frame(time=c(0,0,0), stage=c(1,2,3), density=n0)
for (i in 1:20) {
	n0 <- m%*%n0
	temp <- rbind(temp, c(i, 1, n0[1]))
	temp <- rbind(temp, c(i, 2, n0[2]))
	temp <- rbind(temp, c(i, 3, n0[3]))
}
