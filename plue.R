logistic <- function (p) exp(p)/(1+exp(p))
logit <- function (p) log(p) - log(1-p)
curve(logistic(x),-9,9)

set.seed(42)
# pop iniciail: juv, ad, total
n <- c(10, 15); n.t <- sum(n)
# x obs: maturados, nascidos, sobrev
obs <- c(3, 2, 23)
# NOVA observacao, coerente com as observacoes acima
new.n <- c(round(obs[3]/n.t*(n[1]-obs[1]))+obs[2],round(obs[3]/n.t*n[2])+obs[1]) 
sigma <- obs[3]/n.t
f <- obs[2]/n[2]
gamma <- obs[1]/n[1]
A <- matrix(c(sigma*(1-gamma),f,sigma*gamma,sigma),byrow=T, ncol=2)
tr <- function (A) return(A[1,1]+A[2,2])
A.to.lambda <- function(A) 1/2*(tr(A) + sqrt((tr(A)^2 - 4*det(A))))
lambda <- A.to.lambda(A)
# Agora vamos brincar de Metropolis!
# jumping distribution
Q <- function (x) rnorm (3, as.numeric(x), c(0.02,0.02, 0.02))
# probability distribution
Met.f <- function (x) {
	t <- try(dbinom(obs[3], n.t, as.numeric(x[1]), log=TRUE) +
				  dbinom(obs[2], n[2], as.numeric(x[2]), log=TRUE) +
				  dbinom(obs[1], n[1], as.numeric(x[3]), log=TRUE))
	if (is.nan(t)) return (-Inf)
	return(t)
}
# Ponto inicial:
x = data.frame(sigma=sigma, f=f, gamma=gamma, nLL = -Met.f(c(sigma,f,gamma)))
atual <- as.numeric(x[1,])
#
#Iteration
#
for (i in 1:9999) {
	novo <- Q(atual)
	alfa <- Met.f(novo) - Met.f(atual)
	if (is.nan(alfa)) alfa = -Inf
	if (alfa >= 0 || runif(1,0,1) < exp(alfa)) { #aceite
		x[nrow(x)+1,] <-  c(novo, -Met.f(novo))
		atual <- novo
	} else {
		x[nrow(x)+1,] <- c(atual, -Met.f(atual))
	}
}
dis <- x
getlambda = function (sigma, f, gamma) A.to.lambda (matrix(c(sigma*(1-gamma), f, sigma*gamma, sigma), ncol=2, byrow=TRUE) )
getlambda = Vectorize(getlambda)
dis <- cbind(dis, lambda = getlambda(dis[,1], dis[,2], dis[,3]))

#save(dis, file="R/mini.Rdata")
#load("R/mini.Rdata")
d <- density(dis[,5])

curve(dbinom(obs[3], n.t, x))
curve(dbinom(obs[2], n[2], x), col = 2, add=TRUE)
curve(dbinom(obs[1], n[1], x), col = 3, add=TRUE)

plot(density(dis[,4]))

dis
