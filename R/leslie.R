source("pse.R")
library(sensitivity)
# numero de classes de idade
# Classe 1:	0-19 semanas
# Classe 2: 20-43
# Classe 3: 44-71
# Classe 4: 72-103
# Classe 5: 104+

N <- 350 # N amostras hipercubo

# parametros de fecundidade:
# tomados de Smith & Trout, 1994, tab 1
fertility <- read.table("fertility.txt", header=TRUE)
# Corrigindo os dados pela duracao de cada classe de idade:
fertility$Class2 <- fertility$Class2 / 52 * 24
fertility$Class3 <- fertility$Class3 / 52 * 28
fertility$Class4 <- fertility$Class4 / 52 * 32
fertility$Class5 <- fertility$Class5 / 52 * 36
f2 <- LHSsample (N, name="f2", qprob=qnorm, sum(fertility$Class2), sd(fertility$Class2))
f3 <- LHSsample (N, name="f3", qprob=qnorm, sum(fertility$Class3), sd(fertility$Class3))
f4 <- LHSsample (N, name="f4", qprob=qnorm, sum(fertility$Class4), sd(fertility$Class4))
f5 <- LHSsample (N, name="f5", qprob=qnorm, sum(fertility$Class5), sd(fertility$Class5))

# probabilidades de sobrevivencia
# Vindas da tabela 2, normal
srv <- read.table("survival.txt", header=FALSE)[,1]

p1 <- LHSsample (N, name="p1", qprob=qnorm, prod(srv[1:5]), sd(srv)/2)
p2 <- LHSsample (N, name="p2", qprob=qnorm, prod(srv[6:11]), sd(srv))
p3 <- LHSsample (N, name="p3", qprob=qnorm, prod(srv[12:14]* srv[14]^4), sd(srv))
p4 <- LHSsample (N, name="p4", qprob=qnorm, srv[14]^8, sd(srv))

# LHScorcorr continua pentelho - preciso melhorar esse codigo!!!
vars <- as.data.frame(LHScorcorr(cbind(f2,f3,f4,f5,p1,p2,p3,p4), eps=10e-6))
colnames(vars) <- c("f2","f3","f4","f5","p1","p2","p3","p4")

# Gera a matriz de Leslie
# Calcula o autovalor dominante
innerModel <- function (f2, f3, f4, f5, p1, p2, p3, p4) {
# A matriz media dah um lambda de 1.1536
	L <- matrix(
			c( 0, f2, f3, f4, f5,
			  p1,  0,  0,  0,  0, 
		       0, p2,  0,  0,  0, 
			   0,  0, p3,  0,  0, 
			   0,  0,  0, p4,  0), nrow=5, ncol=5, byrow=TRUE)
	return(as.numeric(eigen(L)$values[2]))
}

# Para realizar a analise FAST, eh necessario escrever um "wrapper"
# que realiza mapply
model <- function (x) {
		return(mapply(innerModel, x[,1], x[,2],x[,3],x[,4],x[,5],x[,6],x[,7],x[,8]))
}

res <- model(vars)

pdf(file="Leslie2.pdf")
par(mfrow=c(2,2))
# Distribuicao dos lambdas
e <- ecdf(res)
plot(e, do.points=FALSE, main="ECDF")
segments(c(1,0), c(0,e(1)), c(1,1),c(e(1),e(1)), lty=2, col="gray70") 

# Partial rank correlation coefficients
plot(pcc(vars, res, nboot=100, rank=TRUE))
abline(h=0, lty=3)

# Morris Screening Method
morr <- morris (model   = model, 
				factors = c("f2", "f3", "f4", "f5", "p1", "p2", "p3", "p4"),
				r       = c(5,10),
				design  = list(type="simplex", scale.factor=1))
plot(morr, xlim=c(0,1))
title("Morris screening")

# Analise eFAST
fast <- fast99 (model   = model, 
				factors = c("f2", "f3", "f4", "f5", "p1", "p2", "p3", "p4"),
				n       = N,
				q       = "qnorm", 
				q.arg   = list(
							   list(mean=sum(fertility$Class2), sd=sd(fertility$Class2)),
							   list(mean=sum(fertility$Class3), sd=sd(fertility$Class3)),
							   list(mean=sum(fertility$Class4), sd=sd(fertility$Class4)),
							   list(mean=sum(fertility$Class5), sd=sd(fertility$Class5)),
							   list(mean=prod(srv[1:5]), sd=sd(srv)/2),
							   list(mean=prod(srv[6:11]), sd=sd(srv)),
							   list(mean=prod(srv[12:14])*srv[14]^4, sd=sd(srv)),
							   list(mean=srv[14]^8, sd=sd(srv))
						      )
				)

plot(fast)
title(main="eFAST")

dev.off()
