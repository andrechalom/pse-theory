source("pse.R")
library(sensitivity)
# Hipercubo gerado de tamanho 650
load("leslie.Rdata")

################# PARTE 1: Independente de Densidade
ivars <- vars[,1:14]

# Gera a matriz de Leslie
# Calcula o autovalor dominante
LeslieIndep <- function (p11, p17, p21, p22, p32, p33, p43, p44, p54, p55, p65, p66, p76, p77) {
	L <- matrix(
			c(p11,   0,   0,   0,   0,   0, p17,
			  p21, p22,   0,   0,   0,   0,   0,
		        0, p32, p33,   0,   0,   0,   0,
			    0,   0, p43, p44,   0,   0,   0,
			    0,   0,   0, p54, p55,   0,   0,
			    0,   0,   0,   0, p65, p66,   0,
			    0,   0,   0,   0,   0, p76, p77), nrow=7, ncol=7, byrow=TRUE)
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
prcc <- pcc(ivars, ires, nboot=100, rank=TRUE)
# Sorted by PRCC:
order(- iprcc$PRCC$original) -> O
iprcc$X <- iprcc$X[,O]
iprcc$PRCC <- iprcc$PRCC[O,]
# Morris Screening Method
imorr <- morris (model   = model, 
				factors = names(ivars),
				r       = 20,
				design  = list(type="simplex", scale.factor=1))
# Analise eFAST
qarg <- list()
for (i in 1:dim(ivars)[2]) 
		qarg[[i]] <- list(mean=mean(ivars[,i]), sd=sd(ivars[,i]))
fast <- fast99 (model   = model, 
				factors = names(ivars),
				n       = 350,
				q       = "qnorm", 
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
LeslieDep <- function (p11, p17, p22, p32, p33, p43, p44, p54, p55, p65, p66, p76, p77, Gm, a, z) {
		assign("Gl", get("Gl",envir=myenv)+1, envir=myenv)
		Np <- c(100,rep(1, 6))
		epsilon = 10e-4
		exit = FALSE
		nIter <- 0 
		while (exit == FALSE) {
				nIter <- nIter + 1
				p21 <- densdep(Gm, a, Np[1], z, Np[7])
				L <- matrix(
							c(p11,   0,   0,   0,   0,   0, p17,
							  p21, p22,   0,   0,   0,   0,   0,
							  0, p32, p33,   0,   0,   0,   0,
							  0,   0, p43, p44,   0,   0,   0,
							  0,   0,   0, p54, p55,   0,   0,
							  0,   0,   0,   0, p65, p66,   0,
							  0,   0,   0,   0,   0, p76, p77), nrow=7, ncol=7, byrow=TRUE)
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
dprcc <- pcc(dvars, dres, nboot=100, rank=TRUE)
# Sorted by PRCC:
order(- dprcc$PRCC$original) -> O
dprcc$X <- dprcc$X[,O]
dprcc$PRCC <- dprcc$PRCC[O,]
# Morris Screening Method
dmorr <- morris (model   = model, 
				factors = names(dvars),
				r       = 20,
				design  = list(type="simplex", scale.factor=1))
# Analise eFAST
qarg <- list()
for (i in 1:dim(dvars)[2]) 
		qarg[[i]] <- list(mean=mean(dvars[,i]), sd=sd(dvars[,i]))
assign("Gl",0, envir=myenv)
dfast <- fast99 (model   = model, 
				factors = names(dvars),
				n       = 66,
				q       = "qnorm", 
				q.arg   = qarg
				)

save(dres, dprcc, dmorr, dfast, dvars, file="Dependent.Rdata")
