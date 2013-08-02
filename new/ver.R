matriz = matrix(c(0, 2, 0.6, 0), ncol=2, byrow=TRUE)
lambda = eigen(matriz)$values[1]
w = eigen(matriz)$vectors[,1]
# left eigenvectors are a bit more tricky:
eL1 = eigen(t(matriz))$vectors
v = eL1[,1]
#Classical sensitivities. M = a21; F = a12
dldM = v[1]*w[2] / v%*%w
elasM = 0.6 / lambda * dldM
dldF = v[2]*w[1] / v%*%w
elasF = 2 / lambda * dldM

#Likelihood-based sensitivities: gerando a distribuicao para M
Mh = function (p) dbinom(6000, 10000, p) / dbinom(6000, 10000, 0.6)

Ph = function (p) 3.65*dpois(2, p)

# Vou gerar N samplesi "tentativos" dessa distribuicao:
n = 1000000
dis = data.frame(M=NA, P=NA)
for (i in 1:n) {
	mysampleM = runif(1,0,1)
	myprobM = runif(1,0,1)
	mysampleP = runif(1,0,30)
	myprobP = runif(1,0,1)
	if (myprobM < Mh(mysampleM) && myprobP < Ph(mysampleP)) dis = rbind(dis, c(mysampleM, mysampleP))
}
dis = dis[-1,]
dim(dis)[1]/n

# para comparar com as dists originais:
plot(density(dis$M))
mult = integrate(Mh,0,1)$value
M = function (p) Mh(p)/mult
curve(M(x), col='red', add=TRUE)

plot(density(dis$P))
mult = integrate(Ph,0,Inf)$value
P = function (p) Ph(p)/mult
curve(P(x), col='red', add=TRUE)

# Correlacao entre as variaveis (nesse caso, espuria)
cor(dis$M, dis$P)

# Analise de incerteza sobre lambda
getlambda = function (m, p) eigen( matrix(c(0, p, m, 0), ncol=2, byrow=TRUE) )$values[1]
getlambda = Vectorize(getlambda)

dis <- cbind(dis, lambda = getlambda(dis$M, dis$P))

plot(ecdf(dis$lambda))

# Analise de sensibilidade sobre lambda
source("pse.R")
corPlot(dis[,c(1,2)], dis[,3])

plot(pcc(dis[,c(1,2)], dis[,3], rank=TRUE))

sl<-estim.slr(dis[,c(3,1,2)])
  f <- apply (dis[,c(1,2)], 2, mean, na.rm=T)
  m <- mean(dis[,3])
  sl * f / m

# Verossimilhanca para gerar as	distribuicoes:
# primeiro censo, s=10, total=30
# segundo censo, s=7, total=25

# binomial com uma observacao, s = 10, total = 30
L = function (p) { return (931395465 * p^10 * (1-p)^20) }
curve(L(x), add=T)
curve (dbinom(10, 30, p=x))
D = function (p) dbinom(10, 30, p)
mult = integrate(D,0,1)$value
cD = function (p) integrate(D, 0, p)$value/mult
cD = Vectorize(cD)
curve(cD)

# veros. cumulativa da binomial:
cL = function (p) { return (integrate(L, 0, p)$value)}
cL = Vectorize(cL)
curve(cL(x))

# Usando dois censos independentes, com as hipoteses FORTES de que
# 1) a probabilidade de transicao se mantem constante a cada ano
# 2) a probabilidade de um evento com um dado individuo eh independente a cada ano

# Segundo censo
L2 = function (p) {return (12498200 * p^7 * (1-p)^18) } 
cL2 = function (p) { return (integrate(L2, 0, p)$value)}
cL2 = Vectorize(cL2)
curve(cL2(x), col="blue", add=T)

# Censos combinados
comboL = function(p) {L(p)*L2(p)/3.045803}
combo = function (p) {return (integrate(comboL,0,p)$value)}
combo = Vectorize(combo)
curve(combo, col="red", add=T)

# Ou podemos aproveitar que se X1 e X2 ~ binom, X1+X2 eh binom tb
somaL = function (p) {return (3821903815930200 * p^17 * (1-p)^(55-17)) }
soma = function (p) {return (integrate(somaL,0,p)$value)}
soma = Vectorize(soma)
curve(soma, col="orange", add=T)
# Vejam, eh igual =D
