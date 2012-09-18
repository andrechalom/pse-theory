source("pse.R")
source("getData.R")
library(sensitivity)
library(msm)
# Hipercubo gerado de tamanho 650
load("leslie.Rdata")

################# PARTE 1: Independente de Densidade
ivars <- vars[,1:14]

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
model <- function (x) {
		return(mapply(LeslieIndep, x[,1], x[,2],x[,3],x[,4],x[,5],x[,6],x[,7],x[,8], x[,9], x[,10], x[,11], x[,12], x[,13], x[,14]))
}
# Roda o modelo e faz as analises:
ires <- model(ivars)

# Partial rank correlation coefficients
iprcc <- pcc(ivars, ires, nboot=1000, rank=TRUE)
# Sorted by PRCC:
order(- iprcc$PRCC$original) -> O
iprcc$X <- iprcc$X[,O]
iprcc$PRCC <- iprcc$PRCC[O,]

# Morris Screening Method
maxes <- apply(ivars, 2, max)
maxes[maxes<1] <- 1
mines <- apply(ivars, 2, min)
mines[mines>0] <- 0
imorr <- morris (model   = model, 
				factors = names(ivars),
				r       = 2,
				design  = list(type="simplex", scale.factor=1),
				binf = mines,
				bsup = maxes
				)

# Analise eFAST

qarg <- list(
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

ifast <- fast99 (model   = model, 
				factors = names(ivars),
				n       = 66,
				q       = "qtnorm",
				q.arg   = qarg
				)

save(ires, iprcc, imorr, ifast, ivars, file="Independent.Rdata")

#################### PARTE 2: DEPENDENTE DE DENSIDADE

dvars <- vars[,c(1,2,4:17)]

# Equacao de dependencia da densidade:
densdep <- function (Gm, a, N1, z, N7) {
		rho <- 25
		G1 <- Gm / (1+a*N1) * exp(-z/rho*N7)
		return (G1)
}
LeslieDep <- function (s1, F7, s2, g2, s3, g3, s4, g4, s5, g5, s6, g6, s7, gm, a, z) {
		assign("Gl", get("Gl",envir=myenv)+1, envir=myenv)
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
	print(c(get("Gl", envir=myenv), nIter))
	return (result);
}
model <- function (x) {
	return(mapply(LeslieDep, x[,1], x[,2],x[,3],x[,4],x[,5],x[,6],x[,7],x[,8], x[,9], x[,10], x[,11], x[,12], x[,13], x[,14], x[,15], x[,16]))
}

myenv <- new.env()
assign("Gl",0, envir=myenv)
dres <- model(dvars)

# Partial rank correlation coefficients
dprcc <- pcc(dvars, dres, nboot=1000, rank=TRUE)
# Sorted by PRCC:
order(- dprcc$PRCC$original) -> O
dprcc$X <- dprcc$X[,O]
dprcc$PRCC <- dprcc$PRCC[O,]
# Morris Screening Method
# Corrige (bug?) com "a" variando de 0 a 0.005
maxes <- apply(dvars, 2, max)
maxes[maxes<1] <- 1
mines <- apply(dvars, 2, min)
mines[mines>0] <- 0
dmorr <- morris (model   = model, 
				factors = names(dvars),
				r       = 2,
				design  = list(type="simplex", scale.factor=1),
				binf = mines,
				bsup = maxes
				)

# Analise eFAST
qarg <- list(
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

assign("Gl",0, envir=myenv); 
dfast <- fast99 (model   = model, factors = names(dvars), n = 66,				
				q       = "qtnorm", 
				q.arg   = qarg
				)

save(dres, dprcc, dmorr, dfast, dvars, file="Dependent.Rdata")
