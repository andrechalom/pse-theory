# Exemplos de utilizacao do pse:
source("pse.R")
# Modelo numerico: crescimento logistico simples
modelRun <- function (Xo, r, K, Time) {
	X <- Xo
	for (i in 0:Time) {
		X <- r*X[length(X)]*(1-X[length(X)]/K)
	}
	return (X)
}
print(modelRun(1,2,3,4))
# Descrevendo o espaco de parametros e gerando o hipercubo
# Total de amostras sera 100:
N <- 100
# r eh uniformemente distribuido entre 0 e 2, semelhante para os demais:
# (Samples para outras distribuicoes serao implementados posteriormente!)
r <- sampleunif(N, 0, 2,"r")
k <- sampleunif(N, 10, 50,"k")
x <- sampleunif(N,1,10,"x")
Time <- sampleunif(N,50,150,"T")
# Rodando o modelo para os diferentes xo, r e k
res <- mapply(modelRun, x, r, k, Time)
X <- data.frame(x,r,k,res)

# Plots e analises
# Plot de sobrevivencia
testPlot(list(r,x,Time,k), test=res>0.0001, mar=c(4,4,1,2))

# Plot da populacao final
Max <- max(X$res)
Min <- min(X$res)
plot(0,0, xlim=limits(r),ylim=limits(k),xlab="Valores de r",ylab="Valores de K",main="Pop Final")
points(X$r,X$k,col=rgb((X$res-Min)/(Max-Min),0,0),pch=19)

# Geramos um modelo linear para prever a populacao em qualquer ponto
# Ele claramente falha ao identificar a superficie de sobrevivencia
modellm <- lm (X$res ~ X$r + X$k)
print (modellm)
print (summary(modellm))
rpred <- rep(1:N/N*(2-0) + 0 - (1/N),N)
kpred <- rep(1:N/N*(50-10) + 10 - 1/N*(50-10)/2, each=N)
prd <- modellm$coefficients[1] + modellm$coefficients[2]*rpred + modellm$coefficients[3]*kpred
# Normalizacao
prd <- (prd-Min) / (Max-Min)
prd[prd < 0] <- 0
prd[prd > Max ] <- Max
plot(0,0, xlim=c(0,2),ylim=c(10,50),xlab="Valores de r",ylab="Valores de K",main="Pop Final Predita E Realizada") 
points(rpred,kpred,col=rgb(prd,0,0),pch=15)
points(X$r,X$k,col=rgb((X$res-Min)/(Max-Min),0,0),pch=23)

# Mais um modelo linear, agora restrito a regiao de sobrevivencia
Xcor <- X[(X$res > 0.001) ,]
modellm <- lm (Xcor$res ~ Xcor$r + Xcor$k)
print (modellm)
print (summary(modellm))
rpred <- rep(1:N/N*(2-1) + 1 - (1/N),N)
kpred <- rep(1:N/N*(50-10) + 10 - 1/N*(50-10)/2, each=N)
prd <- modellm$coefficients[1] + modellm$coefficients[2]*rpred + modellm$coefficients[3]*kpred
# Normalizacao, jogando fora < 0 
prd <- (prd-Min) / (Max-Min)
prd[prd < 0] <- 0
prd[prd > Max ] <- Max
plot(0,0, xlim=c(0,2),ylim=c(10,50),xlab="Valores de r",ylab="Valores de K",main="Pop Final Predita E Realizada") 
points(rpred,kpred,col=rgb(prd,0,0),pch=15)
points(X$r,X$k,col=rgb((X$res-Min)/(Max-Min),0,0),pch=23)

# Correlacao entre k e r
print(cor(k,r))
