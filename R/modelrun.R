### Exemplos de utilizacao do pse:
source("pse.R")
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
N <- 100
# r eh uniformemente distribuido entre 0 e 2
# os outros parametros seguem logica semelhante
r <- LHSsample(N, "r", qunif, 0.25, 2)
k <- LHSsample(N, "k", qunif, 10, 50)
x <- LHSsample(N, "x", qunif, 1, 10)
# O tempo final eh distribuido com forma normal de media = 100 e sd = 20
Time <- LHSsample(N, "T", qnorm, 100, 40)
apply((cbind(r,k,Time,x)), 2, mean)
apply((cbind(r,k,Time,x)), 2, sd)
# Correlacao entre as variaveis geradas. Veja que ela nao eh nec. 0
M <- cor(cbind(r,k,Time,x))
M
max(abs(M[M!=1]))
# Correcao das variaveis para apresentarem correlacao 0
#  newvars <- LHScorcorr (cbind(r,k,Time,x))
# Correcao das variaveis para r e k terem correlacao negativa
MyCor <- matrix(c(   1,-0.1,   -0.2,   -0.3,
		  0.1,   1,   -0.4,   -0.5,
		     0.2,   0.3,   1,   0,
		     0.4,   0.5,   0,   1),4,4)
newvars <- LHScorcorr (cbind(r,k,Time,x), MyCor)
print(MyCor);
print(cor(newvars))

# Uso a notacao [1:N] para manter os atributos - **TODO** refazer pra ficar menos pentelho
k[1:N] <- newvars[,2]
Time[1:N] <- newvars[,3]
x[1:N] <- newvars[,4]
M <- cor(cbind(r,k,Time,x))
M
max(abs(M[M!=1]))

### Rodando o modelo para os diferentes parametros
# O objeto res contem a saida do modelo, no nosso caso, a populacao final
res <- mapply(modelRun, x, r, k, Time)

