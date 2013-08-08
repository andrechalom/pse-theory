library(sensitivity)
source("pse.R")
library(msm) # Para a funcao qtnorm
# Entradas da matriz:
# Vindos de Silva Matos, Freckleton & Watkinson 1999
dados <- read.csv("leslie.csv", header=TRUE)
dados$HECTARE = dados$HECTARE/40
stasis <- dados[dados$TO == dados$FROM & dados$FROM != 7,]
growth <- dados[dados$TO == (dados$FROM +1),]
# Taxa de fecundidade
fertility<- dados$VALUE[dados$TO==1 & dados$FROM==7]
fertilitymean <- mean(dados$VALUE[dados$TO==1 & dados$FROM==7]) 
fertilitysd <- sd(dados$VALUE[dados$TO==1 & dados$FROM==7])
# Calculamos as taxas de sobrevivencia:
sobrev <- data.frame(stasis$FROM, stasis$TO, stasis$YEAR, stasis$VALUE+growth$VALUE, stasis$HECTARE)
colnames(sobrev) <- c("FROM", "TO", "YEAR", "VALUE", "HECTARE")
# Taxa de crescimento real:
realgrowth <- data.frame(stasis$FROM, stasis$TO, stasis$YEAR, (sobrev$VALUE-stasis$VALUE)/sobrev$VALUE, stasis$HECTARE)
colnames(realgrowth) <- c("FROM", "TO", "YEAR", "VALUE", "HECTARE")
realgrowthmeans<-diag(tapply(realgrowth$VALUE, INDEX=list(realgrowth$TO, realgrowth$FROM), FUN=mean))
realgrowthsds<-diag(tapply(realgrowth$VALUE, INDEX=list(realgrowth$TO, realgrowth$FROM), FUN=sd))
# Adicionando de volta a sobrevivencia da classe 7...
sobrev <- rbind(sobrev, dados[dados$TO == 7 & dados$FROM == 7,])
sobrevmeans<-diag(tapply(sobrev$VALUE, INDEX=list(sobrev$TO, sobrev$FROM), FUN=mean))
sobrevsds<-diag(tapply(sobrev$VALUE, INDEX=list(sobrev$TO, sobrev$FROM), FUN=sd))
# N. total de arvores marcadas:
arvores <- tapply(sobrev$HECTARE, INDEX=list(sobrev$TO), FUN=sum)
# Media de gm retirada do paper, sd chutado (impossivel reconstruir)
gmmean <- 0.486/(0.486+mean(dados$VALUE[dados$TO==1 & dados$FROM==1])) 

# Analise grafica dos parametros
plotgrowth <- function() {
plot(realgrowth$VALUE~realgrowth$FROM, pch=10, ylim=c(0, 0.5), cex=sqrt(realgrowth$HECTARE)/3, col=rep(c('red','green', 'blue'), each=6))
arrows(1:6, realgrowthmeans, 1:6, realgrowthmeans+realgrowthsds, length=0.15, angle=90)
arrows(1:6, realgrowthmeans, 1:6, realgrowthmeans-realgrowthsds, length=0.15, angle=90)
}
plotgrowth()

plotsobrev <- function() {
	sobrev<-sobrev[order(sobrev$YEAR),]
plot(sobrev$VALUE~sobrev$FROM, pch=10, ylim=c(0.6, 1.0), cex=sqrt(sobrev$HECTARE)/4, col=rep(c('red','green','blue'), each=7))
arrows(1:7, sobrevmeans, 1:7, sobrevmeans+sobrevsds, length=0.15, angle=90)
arrows(1:7, sobrevmeans, 1:7, sobrevmeans-sobrevsds, length=0.15, angle=90)
}
plotsobrev()

####### Escolha de modelos
library(bbmle)
source("plot-profmle.r")
# Fecundidade: todos os anos igual:
binomNLL<- function(a){
		lambda=exp(a)
	-sum(dpois(fertility,lambda=lambda, log=TRUE))
}
F1 <- mle2(binomNLL, start=list(a=1.34)) # Dah uma probabilidade ~0.79
#plot.profmle(profile(F1))

# Fecundidade com 3 params
binomNLL<- function(a, b, C){
		lambda=exp(a*c(1,0,0)+b*c(0,1,0)+C*c(0,0,1))
	-sum(dpois(fertility,lambda=lambda, log=TRUE))
}
F2 <- mle2(binomNLL, start=list(a=1.34, b=4, C=3)) # Dah uma probabilidade ~0.79
#plot.profmle(profile(F2))

# Fecundidade eh dada por uma poisson
# O parametro da poisson eh dado por uma normal:
binomNLL<- function(mu, sigma){
		lambda=exp(a*c(1,0,0)+b*c(0,1,0)+C*c(0,0,1))
	-sum(dpois(fertility,lambda=lambda, log=TRUE))
}
F2 <- mle2(binomNLL, start=list(a=1.34, b=4, C=3)) # Dah uma probabilidade ~0.79

a1<-coef(F2)[1]
a2<-coef(F2)[2]
a3<-coef(F2)[3]
lambda=exp(a1*c(1,0,0)+a2*c(0,1,0)+a3*c(0,0,1))
-sum(dpois(fertility,lambda=lambda, log=TRUE))

AICtab(F1, F2)

## Modelo 1: TODAS as taxas de sobrevivencia sao iguais
binomNLL<- function(a){
		prob.det=exp(a)/(1 +exp(a))
	-sum(dbinom(round(sobrev$VALUE*sobrev$HECTARE),size=round(sobrev$HECTARE),prob=prob.det, log=TRUE))
}
M1 <- mle2(binomNLL, start=list(a=1.34)) # Dah uma probabilidade ~0.79
#plot.profmle(profile(M1))
p <- exp(coef(M1))/(1+exp(coef(M1)))
plotsobrev()
abline(h = p, lty=3, col='red')

binomNLL<- function(a){
		prob.det=exp(a)/(1 +exp(a))
	-sum(dbinom(round(realgrowth$VALUE*realgrowth$HECTARE),size=round(realgrowth$HECTARE),prob=prob.det, log=TRUE))
}
G1 <- mle2(binomNLL, start=list(a=-0.80)) # Dah uma probabilidade ~0.79
#plot.profmle(profile(G1))
p <- exp(coef(G1))/(1+exp(coef(G1)))
plotgrowth()
abline(h = p, lty=3, col='red')

# Modelo 2: Sobrevivencia cresce com a classe
binomNLL<- function(a, b){
		prob.det=exp(a+b*sobrev$TO)/(1 +exp(a+b*sobrev$TO))
	-sum(dbinom(round(sobrev$VALUE*sobrev$HECTARE),size=round(sobrev$HECTARE),prob=prob.det, log=TRUE))
}
M2 <- mle2(binomNLL, start=list(a=0.62, b=0.62)) 
#plot.profmle(profile(M2))
p1<- coef(M2)[1]
p2<- coef(M2)[2]
plotsobrev()
curve(exp(p1+p2*x)/(1+exp(p1+p2*x)), to=7, add=T, lty=3, col='red')

binomNLL<- function(a, b){
		prob.det=exp(a+b*realgrowth$TO)/(1 +exp(a+b*realgrowth$TO))
	-sum(dbinom(round(realgrowth$VALUE*realgrowth$HECTARE),size=round(realgrowth$HECTARE),prob=prob.det, log=TRUE))
}
G2 <- mle2(binomNLL, start=list(a=0.62, b=0.62)) 
#plot.profmle(profile(G2))
p1<- coef(G2)[1]
p2<- coef(G2)[2]
plotgrowth()
curve(exp(p1+p2*x)/(1+exp(p1+p2*x)), to=7, add=T, lty=3, col='red')

# Modelo 3: saplings tem sobrevivencia diferenciada
binomNLL<- function(a, b){
		p = a+b*(sobrev$TO==1)
		prob.det=exp(p)/(1 +exp(p))
	-sum(dbinom(round(sobrev$VALUE*sobrev$HECTARE),size=round(sobrev$HECTARE),prob=prob.det, log=TRUE))
}
M3 <- mle2(binomNLL, start=list(a=2.21, b=-0.96)) 
#plot.profmle(profile(M3))
p1<- coef(M3)[1]
p2<- coef(M3)[2]
plotsobrev()
segments(c(0.5, 1.5),
		 c(exp(p1+p2)/(1+exp(p1+p2)), exp(p1)/(1+exp(p1))),
		 c(1.5, 7.5),
		   lty=3, col='red')

binomNLL<- function(a, b){
		p = a+b*(realgrowth$TO==1)
		prob.det=exp(p)/(1 +exp(p))
	-sum(dbinom(round(realgrowth$VALUE*realgrowth$HECTARE),size=round(realgrowth$HECTARE),prob=prob.det, log=TRUE))
}
G3 <- mle2(binomNLL, start=list(a=2.21, b=-0.96)) 
#plot.profmle(profile(G3))
p1<- coef(G3)[1]
p2<- coef(G3)[2]
plotgrowth()
segments(c(0.5, 1.5),
		 c(exp(p1+p2)/(1+exp(p1+p2)), exp(p1)/(1+exp(p1))),
		 c(1.5, 7.5),
		   lty=3, col='red')

# Modelo 4: Sobrev diferenciada por ano
binomNLL<- function(a, b, c){
		p = a*(sobrev$YEAR==1991)+b*(sobrev$YEAR==1992)+c*(sobrev$YEAR==1993)
		prob.det=exp(p)/(1 +exp(p))
	-sum(dbinom(round(sobrev$VALUE*sobrev$HECTARE),size=round(sobrev$HECTARE),prob=prob.det, log=TRUE))
}
M4 <- mle2(binomNLL, start=list(a=1.17, b=1, c=1.76))
#plot.profmle(profile(M4))
p1<- coef(M4)[1]
p2<- coef(M4)[2]
p3<- coef(M4)[3]
plotsobrev()
segments(rep(0.5, 3),
		 c(exp(p1)/(1+exp(p1)), exp(p2)/(1+exp(p2)), exp(p3)/(1+exp(p3))),
		 rep(7.5, 3),
		   lty=3, col=c('red', 'green', 'blue'))

binomNLL<- function(a, b, c){
		p = a*(realgrowth$YEAR==1991)+b*(realgrowth$YEAR==1992)+c*(realgrowth$YEAR==1993)
		prob.det=exp(p)/(1 +exp(p))
	-sum(dbinom(round(realgrowth$VALUE*realgrowth$HECTARE),size=round(realgrowth$HECTARE),prob=prob.det, log=TRUE))
}
G4 <- mle2(binomNLL, start=list(a=1.17, b=1, c=1.76))
#plot.profmle(profile(G4))
p1<- coef(G4)[1]
p2<- coef(G4)[2]
p3<- coef(G4)[3]
plotgrowth()
segments(rep(0.5, 3),
		 c(exp(p1)/(1+exp(p1)), exp(p2)/(1+exp(p2)), exp(p3)/(1+exp(p3))),
		 rep(7.5, 3),
		   lty=3, col=c('red', 'green', 'blue'))

# Modelo 5: Sobrev diferenciada por ano E classe
binomNLL<- function(a, b, c, d){
		p = a*(sobrev$YEAR==1991)+b*(sobrev$YEAR==1992)+c*(sobrev$YEAR==1993)+d*sobrev$TO
		prob.det=exp(p)/(1 +exp(p))
	-sum(dbinom(round(sobrev$VALUE*sobrev$HECTARE),size=round(sobrev$HECTARE),prob=prob.det, log=TRUE))
}
M5 <- mle2(binomNLL, start=list(a=0.44, b=0.26, c=1, d=0.64)) # Dah uma probabilidade ~0.79
#plot.profmle(profile(M5))
p1<- coef(M5)[1]
p2<- coef(M5)[2]
p3<- coef(M5)[3]
p4<- coef(M5)[4]
plotsobrev()
curve(exp(p1+p4*x)/(1+exp(p1+p4*x)), to=7, add=T, lty=3, col='red')
curve(exp(p2+p4*x)/(1+exp(p2+p4*x)), to=7, add=T, lty=3, col='green')
curve(exp(p3+p4*x)/(1+exp(p3+p4*x)), to=7, add=T, lty=3, col='blue')
# Salvando para a posteridade
meanp<- mean(c(p1,p2,p3))
sdp<- sd(c(p1,p2,p3))
classp <- p4

binomNLL<- function(a, b, c, d){
		p = a*(realgrowth$YEAR==1991)+b*(realgrowth$YEAR==1992)+c*(realgrowth$YEAR==1993)+d*realgrowth$TO
		prob.det=exp(p)/(1 +exp(p))
	-sum(dbinom(round(realgrowth$VALUE*realgrowth$HECTARE),size=round(realgrowth$HECTARE),prob=prob.det, log=TRUE))
}
G5 <- mle2(binomNLL, start=list(a=0.44, b=0.26, c=1, d=0.64)) # Dah uma probabilidade ~0.79
#plot.profmle(profile(G5))
p1<- coef(G5)[1]
p2<- coef(G5)[2]
p3<- coef(G5)[3]
p4<- coef(G5)[4]
plotgrowth()
curve(exp(p1+p4*x)/(1+exp(p1+p4*x)), to=7, add=T, lty=3, col='red')
curve(exp(p2+p4*x)/(1+exp(p2+p4*x)), to=7, add=T, lty=3, col='green')
curve(exp(p3+p4*x)/(1+exp(p3+p4*x)), to=7, add=T, lty=3, col='blue')

# Modelo I GIVE UP:
binomNLL<- function(a, b, c, d, e, f){
		p = a*(realgrowth$TO==1)+b*(realgrowth$TO==2)+c*(realgrowth$TO==3)+d*(realgrowth$TO==4)+e*(realgrowth$TO==5)+f*(realgrowth$TO==6)
		prob.det=exp(p)/(1 +exp(p))
	-sum(dbinom(round(realgrowth$VALUE*realgrowth$HECTARE),size=round(realgrowth$HECTARE),prob=prob.det, log=TRUE))
}
G6 <- mle2(binomNLL, start=list(a=2.21, b=-0.96, c=1, d=1, e=1, f=1)) 
#plot.profmle(profile(G6))

AICtab(M1, M2, M3, M4, M5, weights=TRUE)
AICtab(G1, G2, G3, G4, G5, G6, weights=TRUE)


################# PARTE 1: Independente de Densidade
factors <- c("s1","F7","g1","s2","g2","s3","g3","s4","g4","s5","g5","s6","g6","s7");
q <- rep("qtnorm", length(factors))
r <- rep("rtnorm", length(factors))
q.arg <- list(
			  list(mean=sobrevmeans[1], sd=sobrevsds[1], lower=0, upper=1),
			  list(mean=fertilitymean, sd=fertilitysd, lower=0, upper=500),
			  list(mean=realgrowthmeans[1], sd=realgrowthsds[1], lower=0, upper=1),
			  list(mean=sobrevmeans[2], sd=sobrevsds[2], lower=0, upper=1),
			  list(mean=realgrowthmeans[2], sd=realgrowthsds[2], lower=0, upper=1),
			  list(mean=sobrevmeans[3], sd=sobrevsds[3], lower=0, upper=1),
			  list(mean=realgrowthmeans[3], sd=realgrowthsds[3], lower=0, upper=1),
			  list(mean=sobrevmeans[4], sd=sobrevsds[4], lower=0, upper=1),
			  list(mean=realgrowthmeans[4], sd=realgrowthsds[4], lower=0, upper=1),
			  list(mean=sobrevmeans[5], sd=sobrevsds[5], lower=0, upper=1),
			  list(mean=realgrowthmeans[5], sd=realgrowthsds[5], lower=0, upper=1),
			  list(mean=sobrevmeans[6], sd=sobrevsds[6], lower=0, upper=1),
			  list(mean=realgrowthmeans[6], sd=realgrowthsds[6], lower=0, upper=1),
			  list(mean=sobrevmeans[7], sd=sobrevsds[7], lower=0, upper=1)
			  )   
# Gera a matriz de Leslie
# Calcula o autovalor dominante
LeslieIndep <- function (s1, F7, g1, s2, g2, s3, g3, s4, g4, s5, g5, s6, g6, s7) {
	L <- matrix(
			c(s1*(1-g1),   0,   0,   0,   0,   0, F7,
			  s1*g1, s2*(1-g2),   0,   0,   0,   0,   0,
		        0, s2*g2, s3*(1-g3),   0,   0,   0,   0,
			    0,   0, s3*g3, s4*(1-g4),   0,   0,   0,
			    0,   0,   0, s4*g4, s5*(1-g5),   0,   0,
			    0,   0,   0,   0, s5*g5, s6*(1-g6),   0,
			    0,   0,   0,   0,   0, s6*g6, s7), nrow=7, ncol=7, byrow=TRUE)
	result <- Mod(eigen(L)$values[1])
	return (result);
}
# Para realizar a analise FAST, eh necessario escrever um "wrapper"
# que realiza mapply
IndepModel <- function (x) {
		return(mapply(LeslieIndep, x[,1], x[,2],x[,3],x[,4],x[,5],x[,6],x[,7],x[,8], x[,9], x[,10], x[,11], x[,12], x[,13], x[,14]))
}

# Roda o modelo e faz as analises:
LHS05 <- LHS(IndepModel, factors, 50, q, q.arg)
LHS1 <- LHS(IndepModel, factors, 100, q, q.arg)
s1<-sbma(LHS05, LHS1)
LHS2 <- LHS(IndepModel, factors, 200, q, q.arg)
s2<-sbma(LHS1, LHS2)
LHS3 <- LHS(IndepModel, factors, 300, q, q.arg)
s3<-sbma(LHS2, LHS3)
LHS4 <- LHS(IndepModel, factors, 400, q, q.arg)
s4<-sbma(LHS3, LHS4)
LHS5 <- LHS(IndepModel, factors, 500, q, q.arg)
s5<-sbma(LHS4, LHS5)

corPlot(LHS1)
plotecdf(LHS1)

# Partial rank correlation coefficients
plotprcc(LHS1)
#iprcc <- pcc(LHS5);order(- abs(iprcc$PRCC$original)) -> O;iprcc$X <- iprcc$X[,O];iprcc$PRCC <- iprcc$PRCC[O,]
#plot(iprcc); abline(h=0, lty=2)

# Analise eFAST
# fast1 <- fast99 (model = IndepModel, factors = factors, n = 1*66, q = q, q.arg = q.arg)
# fast2 <- fast99 (model = IndepModel, factors = factors, n = 2*66, q = q, q.arg = q.arg)
# (fs1 <- sbma(fast1, fast2))
# fast4 <- fast99 (model = IndepModel, factors = factors, n = 4*66, q = q, q.arg = q.arg)
# (fs2 <- sbma(fast2, fast4))
# fast8 <- fast99 (model = IndepModel, factors = factors, n = 8*66, q = q, q.arg = q.arg)
# (fs3 <- sbma(fast4, fast8))
# fast16 <- fast99 (model = IndepModel, factors = factors, n = 16*66, q = q, q.arg = q.arg)
# (fs4 <- sbma(fast8, fast16))
# 
# sobol1 <- mysobol(IndepModel, factors, 200, r, q.arg)
# sobol2 <- mysobol(IndepModel, factors, 2*200, r, q.arg)
# (ss1<-sbma(sobol1, sobol2))
# sobol3 <- mysobol(IndepModel, factors, 4*200, r, q.arg)
# (ss2<-sbma(sobol2, sobol3))
# sobol4 <- mysobol(IndepModel, factors, 8*200, r, q.arg)
# (ss3<-sbma(sobol3, sobol4))
# sobol5 <- mysobol(IndepModel, factors, 4*8*200, r, q.arg)
# (ss4<-sbma(sobol4, sobol5))
# sobol6 <- mysobol(IndepModel, factors, 4*4*8*200, r, q.arg)
# (ss5<-sbma(sobol5, sobol6))
# 

density((get.results(LHS2)[,1])) -> d
# Lambda medio:
(mean(d$x) -> mean.d)
# Lambda mais provavel:
(d$x[which(d$y==max(d$y))] -> prob.d)
(prob.d-mean.d)/prob.d*100

#################### PARTE 2: DEPENDENTE DE DENSIDADE
factors <- c("s1","F7","s2","g2","s3","g3","s4","g4","s5","g5","s6","g6","s7","gm","a","z");
N <- 20
q <- rep("qtnorm", length(factors))
r <- rep("rtnorm", length(factors))
q.arg <- list(
			 list(mean=sobrevmeans[1], sd=sobrevsds[1], lower=0, upper=1),
			 list(mean=fertilitymean, sd=fertilitysd, lower=0, upper=500),
			 list(mean=sobrevmeans[2], sd=sobrevsds[2], lower=0, upper=1),
			 list(mean=realgrowthmeans[2], sd=realgrowthsds[2], lower=0, upper=1),
			 list(mean=sobrevmeans[3], sd=sobrevsds[3], lower=0, upper=1),
			 list(mean=realgrowthmeans[3], sd=realgrowthsds[3], lower=0, upper=1),
			 list(mean=sobrevmeans[4], sd=sobrevsds[4], lower=0, upper=1),
			 list(mean=realgrowthmeans[4], sd=realgrowthsds[4], lower=0, upper=1),
			 list(mean=sobrevmeans[5], sd=sobrevsds[5], lower=0, upper=1),
			 list(mean=realgrowthmeans[5], sd=realgrowthsds[5], lower=0, upper=1),
			 list(mean=sobrevmeans[6], sd=sobrevsds[6], lower=0, upper=1),
			 list(mean=realgrowthmeans[6], sd=realgrowthsds[6], lower=0, upper=1),
			 list(mean=sobrevmeans[7], sd=sobrevsds[7], lower=0, upper=1),
			 list(mean=gmmean, sd=0.06, lower=0, upper=1),
			 list(mean=0.01228, sd=0.01228/3, lower=0, upper=1),
			 list(mean=7, sd=1, lower=0, upper=20)
			 )
# Equacao de dependencia da densidade:
densdep <- function (Gm, a, N1, z, N7) {
		rho <- 25
		G1 <- Gm / (1+a*N1) * exp(-z/rho*N7)
		return (G1)
}
LeslieDep <- function (s1, F7, s2, g2, s3, g3, s4, g4, s5, g5, s6, g6, s7, gm, a, z) {
		Np <- c(100,rep(1, 6))
		epsilon = 10e-4
		exit = FALSE
		nIter <- 0 
		while (exit == FALSE) {
				nIter <- nIter + 1
				g1 <- densdep(gm, a, Np[1], z, Np[7])
	L <- matrix(
			c(s1*(1-g1),   0,   0,   0,   0,   0, F7,
			  s1*g1, s2*(1-g2),   0,   0,   0,   0,   0,
		        0, s2*g2, s3*(1-g3),   0,   0,   0,   0,
			    0,   0, s3*g3, s4*(1-g4),   0,   0,   0,
			    0,   0,   0, s4*g4, s5*(1-g5),   0,   0,
			    0,   0,   0,   0, s5*g5, s6*(1-g6),   0,
			    0,   0,   0,   0,   0, s6*g6, s7), nrow=7, ncol=7, byrow=TRUE)
				newp <- L %*% Np
				if (sqrt (sum(newp-Np)^2) < epsilon | nIter == 500000 ) exit = TRUE
				Np <- newp
		}
		     result <- sum(Np)
	#return (Np);
}
DepModel <- function (x) {
	return(mapply(LeslieDep, x[,1], x[,2],x[,3],x[,4],x[,5],x[,6],x[,7],x[,8], x[,9], x[,10], x[,11], x[,12], x[,13], x[,14], x[,15], x[,16]))
}

dLHS05 <- LHS(DepModel,factors, 50, q, q.arg)
dLHS1 <- LHS(DepModel,factors, 100, q, q.arg)
ds1<-sbma(dLHS05, dLHS1)
dLHS2 <- LHS(DepModel,factors, 200, q, q.arg)
ds2<-sbma(dLHS1, dLHS2)
dLHS3 <- LHS(DepModel,factors, 300, q, q.arg)
ds3<-sbma(dLHS2, dLHS3)
dLHS4 <- LHS(DepModel,factors, 400, q, q.arg)
ds4<-sbma(dLHS3, dLHS4)
dLHS5 <- LHS(DepModel,factors, 500, q, q.arg)
ds5<-sbma(dLHS4, dLHS5)

myLHS <- target.sbma(0.8, DepModel, factors, 100, 100, q, q.arg)

corPlot(dLHS1, index.data=2)
corPlot(dLHS1, index.res=7)

plotecdf(dLHS1)
plotecdf(dLHS1, stack=T, index.res=2:7)

# Partial rank correlation coefficients
plotprcc(dLHS1)

plotelast(dLHS5)
# dprcc <- pcc(dLHS5); order(- abs(dprcc$PRCC$original)) -> O; dprcc$X <- dprcc$X[,O]; dprcc$PRCC <- dprcc$PRCC[O,]
# plot(dprcc); abline(h=0, lty=2)

#### PCA dos resultados
# normalizando os dados pelas medias
# m <- apply(get.results(dLHS1), 2, mean)
# n <- outer(rep(1, nrow(get.results(dLHS1))), m)
# r <- get.results(dLHS1)/n
fit <- princomp(get.results(dLHS1), cor=TRUE)
summary(fit) # print variance accounted for
loadings(fit) # pc loadings
plot(fit,type="lines") # scree plot
fit$scores[,1] # the principal component
p <- pcc(get.data(dLHS1), fit$scores[,1], rank=T)
plot(p)

COMOFAS PRA INTERPRETAR ISSO???
# dfast1 <- fast99 (model = DepModel, factors = factors, n = 1*66, q = q, q.arg = q.arg)
# dfast2 <- fast99 (model = DepModel, factors = factors, n = 2*66, q = q, q.arg = q.arg)
# (dfs1 <- sbma(dfast1, dfast2))
# dfast4 <- fast99 (model = DepModel, factors = factors, n = 4*66, q = q, q.arg = q.arg)
# (dfs2 <- sbma(dfast2, dfast4))
# dfast8 <- fast99 (model = DepModel, factors = factors, n = 8*66, q = q, q.arg = q.arg)
# (dfs3 <- sbma(dfast4, dfast8))
# dfast16 <- fast99 (model = DepModel, factors = factors, n = 16*66, q = q, q.arg = q.arg)
# (dfs4 <- sbma(dfast8, dfast16))
# 
# dsobol <- mysobol(DepModel, factors, 4*8*200, r, q.arg)
# 
# tabprcc <- data.frame(Size=c("50-100", "100-200", "200-300", "300-400", "400-500"), Independent=c(s1,s2,s3,s4,s5), Dependent=c(ds1,ds2,ds3,ds4,ds5))
#tabfast <- data.frame(cbind(c("66-132", "132-264", "264-528", "528-1056"), rbind(fs1, fs2, fs3, fs4), rbind(dfs1, dfs2, dfs3, dfs4)))
# tabfast <- data.frame(cbind(c("66-132", "132-264", "264-528"), rbind(fs1, fs2, fs3), rbind(dfs1, dfs2, dfs3)))
# colnames(tabfast) <-c("Size ($N_s$)", "Indep. $D_i$", "Indep. $D_t$", "Dep. $D_i$", "Dep. $D_t$")
# load("leslie.Rdata")
# save(LHS5, dLHS5, iprcc, dprcc, fast8, dfast8, dfast4, sobol6, dsobol,  tabprcc, tabfast, file="leslie.Rdata")

# Modelando res como funcao das variaveis mais importantes:
# library(bbmle)
# res <- LHS3@res; vars <- LHS3@data;
# 
# M0 <- mle(function (m0=1.2, s0=0.06) -sum(dnorm(res, mean=m0, sd=s0, log=TRUE)))
# M1 <- mle(function (m0=1.2, mf7=0.0005, s0=0.06) -sum(dnorm(res, mean=m0+mf7*vars$F7, sd=s0, log=TRUE)))
# M2 <- mle(function (m0=1.2, mf7=5e-4,mg2=5e-4, s0=0.06) -sum(dnorm(res, mean=m0+mf7*vars$F7+mg2*vars$g2, sd=s0, log=TRUE)))
# M3 <- mle(function (m0=1.2, mf7=5e-4,mg2=5e-4,mg3=5e-4, s0=0.06) -sum(dnorm(res, mean=m0+mf7*vars$F7+mg2*vars$g2+mg3*vars$g3, sd=s0, log=TRUE)))
# M4 <- mle(function (m0=1.2, mf7=5e-4,mg2=5e-4,mg3=5e-4,mg5=5e-4, s0=0.06) -sum(dnorm(res, mean=m0+mf7*vars$F7+mg2*vars$g2+mg3*vars$g3+mg5*vars$g5, sd=s0, log=TRUE)))
# M5 <- mle(function (m0=1.2, mf7=5e-4,mg2=5e-4,mg3=5e-4,mg5=5e-4,ms7=5e-4, s0=0.06) -sum(dnorm(res, mean=m0+mf7*vars$F7+mg2*vars$g2+mg3*vars$g3+mg5*vars$g5+ms7*vars$s7, sd=s0, log=TRUE)))
# M6 <- mle(function (m0=1.2, mf7=5e-4,mg2=5e-4,mg3=5e-4,mg5=5e-4,ms4=5e-4,ms7=5e-4, s0=0.06) -sum(dnorm(res, mean=m0+mf7*vars$F7+mg2*vars$g2+mg3*vars$g3+mg5*vars$g5+ms4*vars$s4+ms7*vars$s7, sd=s0, log=TRUE)))
# 
# AICtab(M0, M1, M2, M3, M4, M5, M6)
# 

# Dependent variables of ñ: ratio between juveniles and adults:
LeslieDep <- function (s1, F7, s2, g2, s3, g3, s4, g4, s5, g5, s6, g6, s7, gm, a, z) {
		Np <- c(100,rep(1, 6))
		epsilon = 10e-4
		exit = FALSE
		nIter <- 0 
		while (exit == FALSE) {
				nIter <- nIter + 1
				g1 <- densdep(gm, a, Np[1], z, Np[7])
	L <- matrix(
			c(s1*(1-g1),   0,   0,   0,   0,   0, F7,
			  s1*g1, s2*(1-g2),   0,   0,   0,   0,   0,
		        0, s2*g2, s3*(1-g3),   0,   0,   0,   0,
			    0,   0, s3*g3, s4*(1-g4),   0,   0,   0,
			    0,   0,   0, s4*g4, s5*(1-g5),   0,   0,
			    0,   0,   0,   0, s5*g5, s6*(1-g6),   0,
			    0,   0,   0,   0,   0, s6*g6, s7), nrow=7, ncol=7, byrow=TRUE)
				newp <- L %*% Np
				if (sqrt (sum(newp-Np)^2) < epsilon | nIter == 500000 ) exit = TRUE
				Np <- newp
		}
	result <- Np[1]/Np[7];
	return (result);
}
DepModel <- function (x) {
	return(mapply(LeslieDep, x[,1], x[,2],x[,3],x[,4],x[,5],x[,6],x[,7],x[,8], x[,9], x[,10], x[,11], x[,12], x[,13], x[,14], x[,15], x[,16]))
}

dLHS05 <- LHS(DepModel,factors, 50, q, q.arg)
dLHS1 <- LHS(DepModel,factors, 100, q, q.arg)
ds1<-sbma(dLHS05, dLHS1)
dLHS2 <- LHS(DepModel,factors, 200, q, q.arg)
ds2<-sbma(dLHS1, dLHS2)
dLHS3 <- LHS(DepModel,factors, 300, q, q.arg)
ds3<-sbma(dLHS2, dLHS3)
dLHS4 <- LHS(DepModel,factors, 400, q, q.arg)
ds4<-sbma(dLHS3, dLHS4)
dLHS5 <- LHS(DepModel,factors, 500, q, q.arg)
ds5<-sbma(dLHS4, dLHS5)

plot(ecdf(dLHS5))
plot(pcc(dLHS5))


#### Partial correlations:
source("caswell.r")
SRC <-src(dLHS5@data, dLHS5@res)
par(mfrow=c(2,1), mar=c(3,3,1,2))
plot(SRC, main="")
abline(h=0, lty=2)
barplot(dNdtheta, names=c("s1","F7","s2","g2","s3","g3","s4","g4","s5","g5","s6","g6","s7","gm","a","z"));
abline(h=0, lty=2)

# TODOOOOOOO
estim.pcc <- function(data, i = 1:nrow(data)) {  
		  d <- data[i, ]
  p <- ncol(d) - 1
    pcc <- numeric(p)
    for (j in 1:p) {
			    Xtildej.lab <- paste(colnames(d)[c(-1, -j-1)], collapse = "+")
	    lm.Y <- lm(formula(paste(colnames(d)[1], "~", Xtildej.lab)), data = d)
		    lm.Xj <- lm(formula(paste(colnames(d)[j+1], "~", Xtildej.lab)), data = d)
		    pcc[j] <- cor(d[1] - fitted(lm.Y), d[j+1] - fitted(lm.Xj))
			  }
	  pcc
}


pcc <- function(X, y, rank = FALSE, nboot = 0, conf = 0.95) {
		  data <- cbind(Y = y, X)

  if (rank) {
		      for (i in 1:ncol(data)) {
					        data[,i] <- rank(data[,i])
      }
    }
    
    if (nboot == 0) {
			    pcc <- data.frame(original = estim.pcc(data))
	    rownames(pcc) <- colnames(X)

