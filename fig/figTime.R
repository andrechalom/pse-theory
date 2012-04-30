### Exemplos de utilizacao do pse:
source("R/pse.R")
### Modelo numerico: crescimento logistico simples
# O modelo tem 4 parametros: r, K, populacao inicial e tempo final
# Vamos estudar qual a influencia de cada um deles
modelRun <- function (Xo, r, K, Time) {
	X <- Xo
	for (i in 0:Time) {
		X <- r*X[length(X)]*(1-X[length(X)]/K)
	}
	return (X)
}

### Descrevendo o espaco de parametros e gerando o hipercubo
# Definicao do numero total de amostras:
N <- 50

PRCC<- matrix(0,20,3)
colnames(PRCC)<-c("r","k","x")
for (Time in (1:20)) {
# r eh uniformemente distribuido entre 0 e 2
# os outros parametros seguem logica semelhante
r <- LHSsample(N, "r", qunif, 0.25, 2)
k <- LHSsample(N, "k", qunif, 10, 50)
x <- LHSsample(N, "x", qunif, 1, 10)
# Correcao das variaveis para apresentarem correlacao 0
newvars <- LHScorcorr (cbind(r,k,x))
# Correcao das variaveis para r e k terem correlacao negativa
MyCor <- matrix(c(   1,-0.5,   0,
		  -0.5,   1,   0,
		     0,   0,   1),3,3)
#newvars <- LHScorcorr (cbind(r,k,Time,x), MyCor)
# Uso a notacao [1:N] para manter os atributos - **TODO** refazer pra ficar menos pentelho
k[1:N] <- newvars[,2]
x[1:N] <- newvars[,3]

### Rodando o modelo para os diferentes parametros
# O objeto res contem a saida do modelo, no nosso caso, a populacao final
res <- mapply(modelRun, x, r, k, Time)

# Rank transformations:
resr <- rank(res); rr <- rank(r); kr <- rank(k); xr <- rank(x)
#PRCC between res & r
PRCC[Time,1]<-
cor(
    residuals(lm(rr~kr+xr)),
    residuals(lm(resr~kr+xr))
)
#PRCC between res & k
PRCC[Time,2] <-
cor(
    residuals(lm(kr~rr+xr)),
    residuals(lm(resr~rr+xr))
)
#PRCC between res & x
PRCC[Time,3] <-
cor(
    residuals(lm(xr~kr+rr)),
    residuals(lm(resr~kr+rr))
)
}

par(cex=1.3)
plot(0, ylim=c(-1,1), xlim=c(0,20),ylab="PRCC",xlab="Time",type="n")
abline(h=0)
abline(h=0.2, lty=3)
abline(h=-0.2, lty=3)
points(PRCC[,1], pch=1, cex=0.8)
points(PRCC[,2], pch=4, cex=0.8)
points(PRCC[,3], pch=6, cex=0.8)
legend("bottomleft",legend=c("r","k",expression(X[0])), pch=c(1,4,6))
