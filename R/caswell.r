### Cálculos de índices de sensibilidade analiticamente, pelo
# método de Caswell (08, 09, 10)


library(MASS)
########### TRIBOLIUM
########## PARAMETROS ######
b = 6.598
cea = 1.155e-2
cel = 1.209e-2
cpa = 4.7e-3
mua = 7.729e-3
mu1 = 2.055e-1
theta <- c(b, cea, cel, cpa, mua, mu1)
ct <- matrix(c(9, 1, 4.5), ncol=3) # metabolic rates
########### POP EM EQUILIBRIO #####
n<- matrix(c(22.6, 18, 385.2), nrow=3) ## equilibrium n
############ MATRIZ A ########
xp <- exp(-cel*n[1]-cea*n[3])
A.n <- matrix(c(0,0, b*xp, 1-mu1, 0, 0, 0, exp(-cpa*n[3]), 1-mua), nrow=3, byrow=T)
############ AUXILIARES
n.1 <- diag(as.vector(1/n))
Mtheta <- diag(as.vector(theta))
I <- diag(3)
K <- kronecker(t(n), I)
############ DERIVADAS
dvecAdnt <- matrix(c(
					 0,0,0, 
					 0,0,0, 
					 0,0,0, 
					 0,0,0, 
					 0,0,0, 
					 0,0,-cpa*exp(-cpa*n[3]), 
					 -b*cel*xp,0, -b*cea*xp,
					 0,0,0, 
					 0,0,0), nrow=9, byrow=T)
dvecAdtheta <- matrix(c(
				     0,0,0,0,0,0,
				     0,0,0,0,0,-1,
				     0,0,0,0,0,0,
				     0,0,0,0,0,0,
				     0,0,0,0,0,0,
				     0,0,0,-n[3]*exp(-cpa*n[3]),0,0,
				     xp,-b*n[3]*xp,-b*n[1]*xp,0,0,0,
				     0,0,0,0,0,0,
				     0,0,0,0,-1,0), nrow=9, byrow=T)
#sensibilidade
dndtheta <- ginv(I-A.n-K%*%dvecAdnt) %*% K %*% dvecAdtheta
#elasticidade
elas <- n.1 %*% dndtheta %*% Mtheta
dNdtheta <- as.numeric(1/(ct %*% n)) * ct %*% dndtheta %*% Mtheta
barplot(dNdtheta, names=c("b", "cea", "cel", "cpa", "mua", "mu1"))

####### E. EDULIS  ###################
####### PARAMETROS ###################
ct <- matrix(1, ncol=7, nrow=1) # sem ponderacao, pop total
dados <- read.csv("leslie.csv", header=TRUE)
stasis <- dados[dados$TO == dados$FROM & dados$FROM != 7,]
growth <- dados[dados$TO == (dados$FROM +1),]
F7 <- mean(dados$VALUE[dados$TO==1 & dados$FROM==7]) 
# Calculamos as taxas de sobrevivencia:
sobrev <- data.frame(stasis$FROM, stasis$TO, stasis$YEAR, stasis$VALUE+growth$VALUE)
colnames(sobrev) <- c("FROM", "TO", "YEAR", "VALUE")
s<-diag(tapply(sobrev$VALUE, INDEX=list(sobrev$TO, sobrev$FROM), FUN=mean))
# A media reportada para p77 eh considerada "alta demais" pelos autores. Correcao arbitraria
s[7] <- mean(dados$VALUE[dados$TO==7 & dados$FROM==7])-0.015 
# Taxa de crescimento real:
realgrowth <- data.frame(stasis$FROM, stasis$TO, stasis$YEAR, (sobrev$VALUE-stasis$VALUE)/sobrev$VALUE)
colnames(realgrowth) <- c("FROM", "TO", "YEAR", "VALUE")
g<-diag(tapply(realgrowth$VALUE, INDEX=list(realgrowth$TO, realgrowth$FROM), FUN=mean))
gm <- 0.486/(0.486+mean(dados$VALUE[dados$TO==1 & dados$FROM==1])) 
a<-0.01228
z <- 7
rho <- 25 ### FIXO!!!
theta <- c(s[1],F7,s[2],g[2],s[3],g[3],s[4],g[4],s[5],g[5],s[6],g[6],s[7],gm,a,z);
############# MODELO #####################
densdep <- function (Gm, a, N1, z, N7) {
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
	return (Np);
}
TheModel <- function (x) {
	return(LeslieDep(x[1], x[2],x[3],x[4],x[5],x[6],x[7],x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15], x[16]))
}
n <- TheModel(theta) ## Populacao em equilibrio
############ MATRIZ A ########
g[1] <- densdep(gm, a, n[1], z, n[7])
A <- matrix(
		c(s[1]*(1-g[1]),   0,   0,   0,   0,   0, F7,
		s[1]*g[1], s[2]*(1-g[2]),   0,   0,   0,   0,   0,
		0, s[2]*g[2], s[3]*(1-g[3]),   0,   0,   0,   0,
		0,   0, s[3]*g[3], s[4]*(1-g[4]),   0,   0,   0,
		0,   0,   0, s[4]*g[4], s[5]*(1-g[5]),   0,   0,
		0,   0,   0,   0, s[5]*g[5], s[6]*(1-g[6]),   0,
		0,   0,   0,   0,   0, s[6]*g[6], s[7]), nrow=7, ncol=7, byrow=TRUE)
############ AUXILIARES
n.1 <- diag(as.vector(1/n))
Mtheta <- diag(as.vector(theta))
I <- diag(7)
K <- kronecker(t(n), I)
############ DERIVADAS
dvecAdnt <- matrix(0, nrow=49, ncol=7) # Apenas as posicoes A11 e A21 dependem de N:
dvecAdnt[1,1] <- s[1]*g[1]*a/(1+a*n[1]) # dA11 / dN1
dvecAdnt[1,7] <- s[1]*g[1]*z/rho
dvecAdnt[2,] <- -dvecAdnt[1,]
dvecAdtheta <- matrix(0, nrow=49, ncol=16)
dvecAdtheta[1,1] <- 1-g[1] # dA11 / ds1
dvecAdtheta[1,14] <- -s[1]*g[1]/gm # dA11  / dgm
dvecAdtheta[1,15] <- s[1]*g[1]*n[1]/(1+a*n[1]) # dA11 / da
dvecAdtheta[1,16] <- s[1]*g[1]*n[7]/rho # dA11 / dz
dvecAdtheta[2,] <- - dvecAdtheta[1,]
dvecAdtheta[2,1] <- g[1]
dvecAdtheta[9,3] <- 1-g[2]; dvecAdtheta[9,4] <- -s[2]; dvecAdtheta[10,3] <- g[2]; dvecAdtheta[10,4] <- s[2]
dvecAdtheta[17,5] <- 1-g[3]; dvecAdtheta[17,6] <- -s[3]; dvecAdtheta[18,5] <- g[3]; dvecAdtheta[10,6] <- s[3]
dvecAdtheta[25,7] <- 1-g[4]; dvecAdtheta[25,8] <- -s[4]; dvecAdtheta[26,7] <- g[4]; dvecAdtheta[26,8] <- s[4]
dvecAdtheta[33,9] <- 1-g[5]; dvecAdtheta[33,10] <- -s[5]; dvecAdtheta[34,9] <- g[5]; dvecAdtheta[34,10] <- s[5]
dvecAdtheta[41,11] <- 1-g[6]; dvecAdtheta[41,12] <- -s[6]; dvecAdtheta[42,11] <- g[6]; dvecAdtheta[42,12] <- s[6]
dvecAdtheta[43,2] <- 1 
dvecAdtheta[49,13] <- 1 
#sensibilidade
dndtheta <- ginv(I-A-K%*%dvecAdnt) %*% K %*% dvecAdtheta
dNdtheta <- ct %*% dndtheta
# barplot(dNdtheta, names=c("s1","F7","s2","g2","s3","g3","s4","g4","s5","g5","s6","g6","s7","gm","a","z"));
#elasticidade
elas <- n.1 %*% dndtheta %*% Mtheta
dNdtheta <- as.numeric(1/(ct %*% n)) * ct %*% dndtheta %*% Mtheta
barplot(dNdtheta, names=c("s1","F7","s2","g2","s3","g3","s4","g4","s5","g5","s6","g6","s7","gm","a","z"));
