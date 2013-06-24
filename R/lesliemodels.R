library(sensitivity)
source("pse.R")
library(msm) # Para a funcao qtnorm
# Entradas da matriz:
# Vindos de Silva Matos, Freckleton & Watkinson 1999
dados <- read.csv("leslie.csv", header=TRUE)
stasis <- dados[dados$TO == dados$FROM & dados$FROM != 7,]
growth <- dados[dados$TO == (dados$FROM +1),]
# Taxa de fecundidade
fertilitymean <- mean(dados$VALUE[dados$TO==1 & dados$FROM==7]) 
fertilitysd <- sd(dados$VALUE[dados$TO==1 & dados$FROM==7])
# Calculamos as taxas de sobrevivencia:
sobrev <- data.frame(stasis$FROM, stasis$TO, stasis$YEAR, stasis$VALUE+growth$VALUE)
colnames(sobrev) <- c("FROM", "TO", "YEAR", "VALUE")
sobrevmeans<-diag(tapply(sobrev$VALUE, INDEX=list(sobrev$TO, sobrev$FROM), FUN=mean))
sobrevsds<-diag(tapply(sobrev$VALUE, INDEX=list(sobrev$TO, sobrev$FROM), FUN=sd))
# A media reportada para p77 eh considerada "alta demais" pelos autores. Correcao arbitraria
sobrevmeans[7] <- mean(dados$VALUE[dados$TO==7 & dados$FROM==7])-0.015 
sobrevsds[7] <- sd(dados$VALUE[dados$TO==7 & dados$FROM==7])
# Taxa de crescimento real:
realgrowth <- data.frame(stasis$FROM, stasis$TO, stasis$YEAR, (sobrev$VALUE-stasis$VALUE)/sobrev$VALUE)
colnames(realgrowth) <- c("FROM", "TO", "YEAR", "VALUE")
realgrowthmeans<-diag(tapply(realgrowth$VALUE, INDEX=list(realgrowth$TO, realgrowth$FROM), FUN=mean))
realgrowthsds<-diag(tapply(realgrowth$VALUE, INDEX=list(realgrowth$TO, realgrowth$FROM), FUN=sd))
# Media de gm retirada do paper, sd chutado (impossivel reconstruir)
gmmean <- 0.486/(0.486+mean(dados$VALUE[dados$TO==1 & dados$FROM==1])) 

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

