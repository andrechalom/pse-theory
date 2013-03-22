### Exemplos de utilizacao do pse:
source("pse.R")
set.seed(42)
### Modelo numerico: crescimento logistico simples
# O modelo tem 4 parametros: r, K, populacao inicial e tempo final
# Vamos estudar qual a influencia de cada um deles
oneRun <- function (r, K, Xo) {
		Time<-4;
	X <- array();
	X[1] <- Xo; #### VAI SER SOBRESCRITO!
	for (i in 1:Time) {
		X[i] <- X[length(X)] + r*X[length(X)]*(1-X[length(X)]/K)
	}
	return (X)
	#     return (X[length(X)])
}
modelRun <- function (dados) {
		mapply(oneRun, dados[,1], dados[,2], dados[,3])
}
### Descrevendo o espaco de parametros e gerando o hipercubo
factors <- c("r", "K", "X0")
res.names <- paste("Time",1:4)
q <- c("qnorm", "qnorm", "qunif")
r <- c("rnorm", "rnorm", "runif")
q.arg <- list( list(mean=1.5, sd=0.5), list(mean=80, sd=1),
				list(min=1, max=50) )
oneRun(1.2, 12, 5)

meuLHS <- LHS(modelRun, factors, 100, q, q.arg, res.names)

matplot(t(get.results(meuLHS)), type='l')
corPlot(meuLHS)
plotecdf(meuLHS)
plotprcc(meuLHS)
