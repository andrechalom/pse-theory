library(sensitivity)
source("pse.R")
library(msm) # Para a funcao qtnorm
# Entradas da matriz:
# Vindos de Costantino et al., 1996
b = 6.598
cea = 1.155e-2
cel = 1.209e-2
cpa = 4.7e-3
mua = 7.729e-3
mu1 = 2.055e-1
factors <- c("b", "cea", "cel", "cpa", "mua", "mu1");
ct <- c(9, 1, 4.5) # metabolic rates
q <- rep("qtnorm", length(factors))
r <- rep("rtnorm", length(factors))
q.arg <- list(
			  list(mean=b, sd=3e-4, lower=0, upper=100),
			  list(mean=cea, sd=3e-4, lower=0, upper=1),
			  list(mean=cel, sd=3e-4, lower=0, upper=1),
			  list(mean=cpa, sd=3e-4, lower=0, upper=1),
			  list(mean=mua, sd=3e-4, lower=0, upper=1),
			  list(mean=mu1, sd=3e-4, lower=0, upper=1))
Tribolium <- function(b, cea, cel, cpa, mua, mu1) {
		Np <- c(20, 20, 300) #Near steady-state
		epsilon = 10e-5
		exit = FALSE
		nIter <- 0 
		while (exit == FALSE) {
				nIter <- nIter + 1
				L <- matrix(c(0, 0, b*exp(-cel*Np[1] - cea*Np[3]),
							1-mu1, 0, 0,
							0, exp(-cpa*Np[3]), 1-mua), nrow=3, byrow=TRUE)
				newp <- L %*% Np
				if (sqrt (sum(newp-Np)^2) < epsilon | nIter == 500000 ) exit = TRUE
				Np <- newp
		}
	result <- ct %*% Np;
	return (result);
}
# Cheking the output, should be 1950:
Tribolium(b, cea, cel, cpa, mua, mu1)
TheModel <- function (x) {
	return(mapply(Tribolium, x[,1], x[,2],x[,3],x[,4],x[,5],x[,6]))
}

LHS2 <- LHS(TheModel,factors, 200, q, q.arg)
LHS4 <- LHS(TheModel,factors, 400, q, q.arg)
s1 <- sbma(LHS2, LHS4)

plot(ecdf(LHS4))

# Partial rank correlation coefficients
prcc <- pcc(LHS4);
plot(prcc); abline(h=0, lty=2)

sobol <- mysobol(TheModel, factors, 4*8*200*8, r, q.arg)
plot(sobol)
save(LHS4, prcc, sobol, file="Tribolium.rdata")
