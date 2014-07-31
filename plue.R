valid <- function(x) (!is.infinite(x)) & (!is.nan(x)) & (!is.na(x)) 

MCMC.internal <- function(LL, start, opts) {
	# jumping distribution
	if (is.null(opts$Q)) {
		Q <- function (x) rnorm (length(factors), as.numeric(x), rep(0.02, length(factors)))
	} else {
		Q <- opts$Q
	}
	current <- as.numeric(start)
	# L will hold our input data
	L = as.data.frame(matrix(nrow = N, ncol = length(factors)))
	L[1,] <- current
	# nLL (soon to be a vector) will hold the -LL of each data point
	nLL = -LL(current)
	#
	#Metropolis/Monte Carlo Iteration
	#
	for (i in 2:N) {
		new <- Q(current)
		alpha <- LL(new) - LL(current)
		if (is.nan(alpha)) alpha = -Inf
		if (alpha >= 0 || runif(1,0,1) < exp(alpha)) { # point accepted
			L[i,] <-  new
			nLL[i] <- - LL(new)
			current <- new
		} else { # we continue on the same point
			L[i,] <- current
			nLL[i] <- -LL(new)
		}
	}
	return (list(L=L, nLL=nLL))
}



# Use PLUE to perform a Profile Likelihood Uncertainty Analysis
# Comments on the parameters below
PLUE <- function(model, factors, N, LL, start, method=c("internal", "mcmc"), opts = list()) {
	# Input validation for common errors and "default" value handling:
	method = match.arg(method)
	my.opts = list(Q=NULL)
	my.opts[names(opts)] <- opts

	if (is.numeric(factors) && length(factors) == 1) {
		factors = paste("I", 1:factors, sep = "")
	} else if (!is.character(factors)) {
		stop("Error in function PLUE: factors should be either a single number or a character vector")
	}
	if (!is.numeric(N) || length(N) != 1) {
		stop("Error in function PLUE: N should be a single number")
	}
	if (class(LL) != "function") {
		stop("Error in function PLUE: nLL should be a function")
	}

	if (method=="internal") {
		data = MCMC.internal(LL, start, opts)
	} else {
		require(mcmc)
		stop("TODO!")
	}
	# Apply the model to the input data
	res <- model(data$L)

	# Now we profile the results to find out 
	# "what are the lower/upper limits to the 0.1-Delta nLL region?"
	# "what are the lower/upper limits to the 0.2-Delta nLL region?"
	# and so on and so on
	
	mmin <- min(data$nLL[valid(data$nLL) ])
	mmax <- max(data$nLL[valid(data$nLL)])
	prof <- seq(mmin, mmax, length.out=100)
	lower = c(); upper = c();
	for (i in 1:100)
	{
		search <- res[data$nLL <= prof[i]]
		search <- search[valid(search)]
		lower[i] <- min(search); upper[i] <- max(search)
	}
	# finally, we subtract the minimum LL to normalize the profile to 0
	prof = prof - mmin

	X <- list(call = match.call(), N = N, model = model, factors = factors, LL = LL, nLL = data$nLL, start = start,
			  opts = opts, data = data$L, res = res, prof = prof, lower = lower, upper = upper)
	class(X) <- "PLUE"
	return(X)
}
	
print.PLUE <- function(plue) {
	stop("TODO!")
}

# Something like a plot.profmle for PLUE
plot.PLUE <- function(plue) {
	plot(0, type='n', xlim =c(min(plue$lower), max(plue$upper)), 
		 ylim=c(0, max(plue$prof)), main="PLUE", xlab="Result", ylab="Delta Likelihood")
	lines(plue$prof~plue$lower)
	lines(plue$prof~plue$upper)
	abline(h=2, lty=3)
}

####################### END PLUE CODE 
# An example of how to use it (see pse.pdf sec. 3.7)
# Helper functions for the model
tr <- function (A) return(A[1,1]+A[2,2])
A.to.lambda <- function(A) 1/2*(tr(A) + sqrt((tr(A)^2 - 4*det(A))))
getlambda = function (sigma, f, gamma) {
	A.to.lambda (matrix(c(sigma*(1-gamma), f, sigma*gamma, sigma), ncol=2, byrow=TRUE) )
}
getlambda = Vectorize(getlambda)
model <- function(x) getlambda(x[,1], x[,2], x[,3])

factors = c("sigma", "f", "gamma")

N = 20000

# pop iniciail: juv, ad, total
n <- c(10, 15); n.t <- sum(n)
# x obs: maturados, nascidos, sobrev.totais
obs <- c(3, 2, 23)
# melhor chute para os parametros
sigma <- obs[3]/n.t
f <- obs[2]/n[2]
gamma <- obs[1]/n[1]
start = c(sigma, f, gamma)
# probability distribution. It's the POSITIVE LL, because I inverted it somewhere.
LL <- function (x) {
	# try... not working properly (still displays a lot of warning messages!)
	t <- try(dbinom(obs[3], n.t, as.numeric(x[1]), log=TRUE) +
				  dbinom(obs[2], n[2], as.numeric(x[2]), log=TRUE) +
				  dbinom(obs[1], n[1], as.numeric(x[3]), log=TRUE))
	return(t)
}

# On a 3 GHz Pentium, N = 10000 takes 7 secs, N = 100000 takes 700, so growth seems to be quadratic
system.time(plue <- PLUE(model, factors, N, LL, start))
plot(plue)
