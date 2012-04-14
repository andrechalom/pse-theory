# Modelo numerico
modelRun <- function (Xo, r, K, Time) {
	X <- Xo
	for (i in 0:Time) {
		X <- rbind (X, r*X[length(X)]*(1-X[length(X)]/K))
	}
	return (X)
}
print(modelRun(1,2,3,4))

# Descrevendo o espaco de parametros e gerando o hipercubo
N <- 100
# r eh uniformemente distribuido entre 0 e 2
rlimits <- 0:N / N*(2-0) + 0
rpos <- array(0)
for (i in 1:N) rpos[i] <- runif(1, rlimits[i], rlimits[i+1])
rpos <- sample(rpos)
# K eh uniformemente distribuido entre 10 e 50
klimits <- 0:N/N*(50-10) + 10
kpos <- array(0)
for (i in 1:N) kpos[i] <- runif(1, klimits[i], klimits[i+1])
kpos <- sample (kpos)

# Rodando o modelo para os diferentes r e k
# Soh me interessa o resultado final
res <- array(0)
for (i in 1:N) res[i] <- modelRun(5, rpos[i], kpos[i], 100)[100]
print(res)
X <- data.frame(rpos,kpos, res)

# Plots e analises
# Plot de sobrevivencia
plot(0,0, xlim=c(0,2),ylim=c(10,50),xlab="Valores de r",ylab="Valores de K",main="Sobrevivencia")
points(X$rpos[X$res>0.001],X$kpos[X$res>0.001], pch=0)
points(X$rpos[X$res<=0.001],X$kpos[X$res<=0.001], pch=1)
segments(1,10,1,50,lty=2)

# Plot da populacao final
Max <- max(X$res)
Min <- min(X$res)
plot(0,0, xlim=c(0,2),ylim=c(10,50),xlab="Valores de r",ylab="Valores de K",main="Pop Final")
points(X$rpos,X$kpos,col=rgb((X$res-Min)/(Max-Min),0,0),pch=19)

# Geramos um modelo linear para prever a populacao em qualquer ponto
# Ele claramente falha ao identificar a superficie de sobrevivencia
modellm <- lm (X$res ~ X$rpos + X$kpos)
print (modellm)
print (summary(modellm))
rpred <- rep(1:N/N*(2-0) + 0 - (1/N),N)
kpred <- rep(1:N/N*(50-10) + 10 - 1/N*(50-10)/2, each=N)
prd <- modellm$coefficients[1] + modellm$coefficients[2]*rpred + modellm$coefficients[3]*kpred
# Normalizacao, jogando fora < 0 
prd <- prd / max(prd)
prd[prd < 0] <- 0
plot(0,0, xlim=c(0,2),ylim=c(10,50),xlab="Valores de r",ylab="Valores de K",main="Pop Final Predita") 
points(rpred,kpred,col=rgb(prd,0,0),pch=19)

# Correlacao entre k e r
print(cor(kpos,rpos))
