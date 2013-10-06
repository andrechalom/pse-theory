# General "non-reproducing juveniles and reproductive adults" model
# A <- | sigma1 (1-gamma)    f      |
#      | sigma1 gamma        sigma2 |
# ASSUMING survival is independent of stage:
# A <- | sigma (1-gamma)    f     | = | a b |
#      | sigma gamma        sigma |   | c d |
# For large animals with one offspring per cycle, 
# f is ~ a proportion of adults that generate offspring
# sigma ~ binom(\theta_1, n1 = n_2+n_3); f <- binom(\theta_2, n_2); gamma <- binom(\theta_3, n_3)
# IF all variables are independent, and both values for n are fixed, given a single sample {x_1, x_2, x_3}
# \mathcal{L} ( \mathbf{\theta} | \mathbf{x} ) = \sum_i log ( {n_i \choose x_i} \theta_^{x_i} (1-\theta_i) ^{n_i-x_i} 
#
# Por outro lado, o resultado do modelo é \lambda, o maior autovalor de A
# obedecendo det | \lambda - a -b          | = \lambda^2 - \lambda tr(A) + det(A) = 0
#                | -c          \lambda - d |
# \lambda = 1/2 * (tr(A) + \sqrt(tr^2(A) - 4 det(A)) )
# queremos saber, dada uma amostra x_i, quais são as sensibilidades de \lambda aos parâmetros \theta_i

# Exemplo usando Metropolis
(N <- c(20, 10, 10))
(obs <- c(18, 3, 1))
(sigma <- obs[1]/N[1])
(f <- obs[2]/N[2])
(gamma <- obs[3]/N[3])
(A <- matrix(c(sigma*(1-gamma),f,sigma*gamma,sigma),byrow=T, ncol=2))
tr <- function (A) return(A[1,1]+A[2,2])
#eigen(A)
(lambda <- 1/2*(tr(A) + sqrt((tr(A)^2 - 4*det(A)))))
#Simulando a pop
n <- c(0,200)
adults <- numeric()
seeds <- numeric()
for (i in 1:10) {
	 n<- A%*%n
	 adults[i] <- n[2]
	 seeds[i] <- n[1]
}
par(mfrow=c(2,1))
plot(seeds)
plot(adults)
# Agora vamos brincar de Metropolis!
# Ponto inicial:
x = data.frame(sigma=0.9, f=0.3, gamma=0.1)
atual <- as.numeric(x[1,])
# jumping distribution
Q <- function (x) rnorm (3, as.numeric(x), c(0.02,0.02, 0.02))
# probability distribution
f <- function (x) dbinom(obs[1], N[1], as.numeric(x[1])) *
				  dbinom(obs[2], N[2], as.numeric(x[2])) *
				  dbinom(obs[3], N[3], as.numeric(x[3]))
#
#Iteration
#
for (i in 1:100000) {
	novo <- Q(atual)
	alfa <- f(novo)/f(atual)
	if (is.nan(alfa)) alfa = 0
	if (alfa >= 1 || runif(1,0,1) < alfa) { #aceite
		x[nrow(x)+1,] <-  novo
		atual <- novo
	} else {
		x[nrow(x)+1,] <- atual
	}
}
dis <- x
# Analise de incerteza sobre lambda
getlambda = function (sigma, f, gamma) eigen( matrix(c(sigma*(1-gamma), f, sigma*gamma, sigma), ncol=2, byrow=TRUE) )$values[1]
getlambda = Vectorize(getlambda)
dis <- cbind(dis, lambda = getlambda(dis$sigma, dis$f, dis$gamma))
# instabilidades numericas...
dis$lambda[dis$lambda<0] <- - dis$lambda[dis$lambda<0]
plot(density(dis$lambda))

# Analise de sensibilidade sobre lambda
source("pse.R")
corPlot(dis[,c(1,2,3)], dis[,4])

plot(pcc(dis[,c(1,2,3)], dis[,4], rank=TRUE))

sl<-estim.slr(dis[,c(4,1,2,3)])
f <- apply (dis[,c(1,2,3)], 2, mean, na.rm=T)
m <- mean(dis[,4])
sl * f / m




############## 
# exemplo babaca de "particao de incerteza"
# y = x1 + x2; presumimos x1 ~ Pois(lambda1); x2 ~ Pois(lambda2), mas não sabemos valores de lambda
# As perguntas relevantes são:
# Incerteza: qual é o valor mais verossímil de y? quanta confiança temos nesse valor?
# Sensibilidade: qual é a contribuição relativa de x1 e x2 para a incerteza de y?
