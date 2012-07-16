########################CORCORR

# Uso: Para forcar correlacao = 0
# newvars <- LHScorcorr(vars)
# vars deve ser uma matrix ou um objeto que pode ser coerced para matrix
# Para aproximar de uma matriz de correlacao dada:
# newvars <- LHScorcorr(vars,COR)

system("R CMD SHLIB corcorr.c")

# NAO MEXA nos parametros l e it
LHScorcorr <- function (vars, COR = matrix(0,dim(vars)[2],dim(vars)[2]), l = 2, eps = 0.0005, it = 1, echo=FALSE, maxIt = 2*sqrt(dim(vars)[1])) {
	N <- dim(vars)[1]; M <- dim(vars)[2]
	# Condicao de parada: terminamos a correcao
	if (l == M + 1) {
		return (vars);
	}
	# Condicao de skip: correlacoes pequenas o suficiente:
	if (max(abs(cor(vars)[l,1:(l-1)] - COR[l,1:(l-1)])) < eps) {
		return (LHScorcorr (vars, COR, l = l + 1, eps = eps, it = 1, echo=echo));
	}
	# Condicao de parada: iteracao maxima atingida para a mesma variavel
	if (it > maxIt) { 
		if (echo==T) print("WARNING: correlacao nao converge para o esperado apos numero maximo de iteracoes");
		return (LHScorcorr (vars, COR, l = l + 1, eps = eps, it = 1, echo=echo));
	}
	if (echo==T) print(paste("INFO: Realizando correcao de correlacao para l =",l,"/",M))
	# Aqui comeca o trabalho para corrigir as cors da var[,l]
	# Chamada externa para procedure em C, para corrigir correlacoes em var[,l]:
	dyn.load("corcorr.so")
	V <- .C("corcorr", vars=as.double(as.matrix(vars)),cor=as.double(COR), N=as.integer(N), M=as.integer(M), l=as.integer(l), FLAGSTOP=as.integer(0))
	vars <- matrix(V$vars, nrow=N, ncol=M)
	if (V$FLAGSTOP == 1) { # O procedimento convergiu
		return(LHScorcorr (vars, COR, l = l + 1, eps = eps, it = 1, echo=echo));
	} else {
		# Variavel corrigida, repete o procedimento com a mesma variavel
		# ateh a correlacao estar pequena o suficiente
		LHScorcorr(vars,COR,l=l,eps=eps,it=it+1, echo=echo)
	}
}
###################################

### Geradores das amostragens e funcoes acessorias:
LHSsample <- function (N, name="", qprob=qunif, ...) {
	pos <- qprob((1:N)/N - 1/N/2, ...)
	pos <- sample(pos);
	attr(pos,"distribution") <- qprob
	attr(pos,"max") <- qprob(0.99, ...)
	attr(pos,"min") <- qprob(0.01, ...)
	attr(pos,"name") <- name
	return (pos);
}

LHSextend <- function (old, ...) {
	qprob <- attr(old, "distribution")
	newpos <- qprob(c((1:N)/N - 5/N/6, (1:N)/N - 1/N/6), ...)
	newpos <- sample(newpos)
	attr(newpos,"distribution") <- qprob
	attr(newpos,"max") <- attr(old, "name")
	attr(newpos,"min") <- attr(old, "name")
	attr(newpos,"name") <- attr(old, "name")
	return (newpos);
}

# Dada uma amostragem gerada por sample*, devolve os limites de plotagem para a variavel
limits <- function (pos) {
	return (c(attr(pos,"min"),attr(pos,"max")))
}

### Plots e analises
sig <- function (form) {
	f <- anova(lm(form))$"Pr(>F)"[1]
	if (f < 0.001) return ("***");
	if (f < 0.01) return ("**");
	if (f < 0.05) return ("*");
	if (f < 0.1) return (".");
	return (" ");
}
oneCorPlot <- function(res, var) {
	plot(var,res)
	abline((lm(res~var)))
	mtext(paste("Cor:",format(cor(var,res), digits=2), sig(res~var)))
}
corPlot <- function (vars, res) {
	nl <- floor(sqrt(length(vars)))
	nc <- ceiling(length(vars)/nl)
	par(mfrow=c(nl,nc))
	for (i in 1:nc) for (j in 1:nl) {
		if ((i-1)*nc + j <= length(vars)) oneCorPlot(res,vars[[(i-1)*nc+j]])
	}
}

# Funcao acessoria - nao eh para uso externo
oneTestPlot <- function(p1, p2, test, convex=FALSE) {
	if (convex) {
		require(spatstat)
		W <- owin(c(min(p1[test]),max(p1[test])), c(min(p2[test]),max(p2[test])))
		P <- ppp(p1[test],p2[test],window=W)
		cTrue <-convexhull(P)
		W <- owin(c(min(p1[!test]),max(p1[!test])), c(min(p2[!test]),max(p2[!test])))
		P <- ppp(p1[!test],p2[!test],window=W)
		cFalse <- convexhull(P)
	}
	plot(0,0, xlim=limits(p1),ylim=limits(p2),xlab=paste("Valores de ",attr(p1,"name")),ylab=paste("Valores de ",attr(p2,"name")))
	if (convex) {
		plot(cTrue, add=T, col=cm.colors(1), density=10)
		plot(cFalse, add=T, col=cm.colors(2)[2], density=10)
	}
	points(p1[test],p2[test], pch='+')
	points(p1[!test],p2[!test], pch='-')
}
# Plota o resultado de um teste logico pelo hipercubo
testPlot <- function(vars, test, convex=FALSE, ...) {
	l <- length(vars)-1
	opar <- par(mfrow=c(l,l), ...)
	for (i in 1:l) {
		for (j in (l+1):2) {
			if (i >= j) {
				plot.new()
			} else {
				oneTestPlot(vars[[j]],vars[[i]],test, convex)
			}
		}
	}
	par(opar)
}

oneGradPlot <- function (p1, p2, res) {
	Max <- max(res)
	Min <- min(res)
	plot(0,0, xlim=limits(p1),ylim=limits(p2),xlab=paste("Valores de ",attr(p1,"name")),ylab=paste("Valores de ",attr(p2,"name")))
	points(p1,p2,col=rgb((res-Min)/(Max-Min),0,0),pch=19)
}

# Plota o resultado da simulacao pelo hipercubo
gradPlot <- function(vars, res, ...) {
	l <- length(vars)-1
	opar <- par(mfrow=c(l,l), ...)
	for (i in 1:l) {
		for (j in (l+1):2) {
			if (i >= j) {
				plot.new()
			} else {
				oneGradPlot(vars[[j]],vars[[i]],res)
			}
		}
	}
	par(opar)
}

onePredPlot <- function(p1, p2, c0, c1, c2, res) {
	p1pred <- rep(1:N/N*(limits(p1)[2]-limits(p1)[1]) + limits(p1)[1] - 1/N*(limits(p1)[2]-limits(p1)[1])/2, N)
	p2pred <- rep(1:N/N*(limits(p2)[2]-limits(p2)[1]) + limits(p2)[1] - 1/N*(limits(p2)[2]-limits(p2)[1])/2, each=N)
	prd <- c0 + p1pred*c1 + p2pred*c2
	# Normalizacao
	Max <- max(res)
	Min <- min(res)
	prd <- (prd-Min) / (Max-Min)
	prd[prd < 0] <- 0
	prd[prd > 1 ] <- 1
	plot(0,0, xlim=limits(p1),ylim=limits(p2),xlab=paste("Valores de ",attr(p1,"name")),ylab=paste("Valores de ",attr(p2,"name"))) 
	points(p1pred,p2pred,col=rgb(prd,0,0),pch=15)
	points(p1,p2,col=rgb((res-Min)/(Max-Min),0,0),pch=23)
}

predPlot <- function(vars, res, model, ...) {
	l <- length(vars)-1
	opar <- par(mfrow=c(l,l), ...)
	model$coefficients -> cf
	for (i in 1:l) {
		for (j in (l+1):2) {
			if (i >= j) {
				plot.new()
			} else {
				rtmp <- res
				for (k in 0:l+1) { if(k != i & k != j) {
					rtmp <- rtmp - cf[k+1] * mean(vars[[k]])  } }
onePredPlot(vars[[j]],vars[[i]],cf[1], cf[j+1],cf[i+1],rtmp)
			}
		}
	}
	par(opar)
}

# Metodo de calcular a PRCC:
# Downloaded from http://www.yilab.gatech.edu/pcor.html
pcor.test <- function(x,y,z,use="mat",method="p",na.rm=T){
	# The partial correlation coefficient between x and y given z
	#
	# pcor.test is free and comes with ABSOLUTELY NO WARRANTY.
	#
	# x and y should be vectors
	#
	# z can be either a vector or a matrix
	#
	# use: There are two methods to calculate the partial correlation coefficient.
	#	 One is by using variance-covariance matrix ("mat") and the other is by using recursive formula ("rec").
	#	 Default is "mat".
	#
	# method: There are three ways to calculate the correlation coefficient, 
	#	    which are Pearson's ("p"), Spearman's ("s"), and Kendall's ("k") methods.
	# 	    The last two methods which are Spearman's and Kendall's coefficient are based on the non-parametric analysis.
	#	    Default is "p".
	#
	# na.rm: If na.rm is T, then all the missing samples are deleted from the whole dataset, which is (x,y,z).
	#        If not, the missing samples will be removed just when the correlation coefficient is calculated.
	#	   However, the number of samples for the p-value is the number of samples after removing 
	#	   all the missing samples from the whole dataset.
	#	   Default is "T".

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z)

	if(use == "mat"){
		p.use <- "Var-Cov matrix"
		pcor = pcor.mat(x,y,z,method=method,na.rm=na.rm)
	}else if(use == "rec"){
		p.use <- "Recursive formula"
		pcor = pcor.rec(x,y,z,method=method,na.rm=na.rm)
	}else{
		stop("\'use\' should be either \"rec\" or \"mat\"!\n")
	}

	# print the method
	if(gregexpr("p",method)[[1]][1] == 1){
		p.method <- "Pearson"
	}else if(gregexpr("s",method)[[1]][1] == 1){
		p.method <- "Spearman"
	}else if(gregexpr("k",method)[[1]][1] == 1){
		p.method <- "Kendall"
	}else{
		stop("\'method\' should be \"pearson\" or \"spearman\" or \"kendall\"!\n")
	}

	# sample number
	n <- dim(na.omit(data.frame(x,y,z)))[1]
	
	# given variables' number
	gn <- dim(z)[2]

	# p-value
	if(p.method == "Kendall"){
		statistic <- pcor/sqrt(2*(2*(n-gn)+5)/(9*(n-gn)*(n-1-gn)))
		p.value <- 2*pnorm(-abs(statistic))

	}else{
		statistic <- pcor*sqrt((n-2-gn)/(1-pcor^2))
  		p.value <- 2*pnorm(-abs(statistic))
	}

	data.frame(estimate=pcor,p.value=p.value,statistic=statistic,n=n,gn=gn,Method=p.method,Use=p.use)
}			

# By using var-cov matrix
pcor.mat <- function(x,y,z,method="p",na.rm=T){

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z)

	if(dim(z)[2] == 0){
		stop("There should be given data\n")
	}

	data <- data.frame(x,y,z)

	if(na.rm == T){
		data = na.omit(data)
	}

	xdata <- na.omit(data.frame(data[,c(1,2)]))
	Sxx <- cov(xdata,xdata,m=method)

	xzdata <- na.omit(data)
	xdata <- data.frame(xzdata[,c(1,2)])
	zdata <- data.frame(xzdata[,-c(1,2)])
	Sxz <- cov(xdata,zdata,m=method)

	zdata <- na.omit(data.frame(data[,-c(1,2)]))
	Szz <- cov(zdata,zdata,m=method)

	# is Szz positive definite?
	zz.ev <- eigen(Szz)$values
	if(min(zz.ev)[1]<0){
		stop("\'Szz\' is not positive definite!\n")
	}

	# partial correlation
	Sxx.z <- Sxx - Sxz %*% solve(Szz) %*% t(Sxz)
	
	rxx.z <- cov2cor(Sxx.z)[1,2]

	rxx.z
}

# By using recursive formula
pcor.rec <- function(x,y,z,method="p",na.rm=T){
	# 

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z)

	if(dim(z)[2] == 0){
		stop("There should be given data\n")
	}

	data <- data.frame(x,y,z)

	if(na.rm == T){
		data = na.omit(data)
	}

	# recursive formula
	if(dim(z)[2] == 1){
		tdata <- na.omit(data.frame(data[,1],data[,2]))
		rxy <- cor(tdata[,1],tdata[,2],m=method)

		tdata <- na.omit(data.frame(data[,1],data[,-c(1,2)]))
		rxz <- cor(tdata[,1],tdata[,2],m=method)

		tdata <- na.omit(data.frame(data[,2],data[,-c(1,2)]))
		ryz <- cor(tdata[,1],tdata[,2],m=method)

		rxy.z <- (rxy - rxz*ryz)/( sqrt(1-rxz^2)*sqrt(1-ryz^2) )
		
		return(rxy.z)
	}else{
		x <- c(data[,1])
		y <- c(data[,2])
		z0 <- c(data[,3])
		zc <- as.data.frame(data[,-c(1,2,3)])

		rxy.zc <- pcor.rec(x,y,zc,method=method,na.rm=na.rm)
		rxz0.zc <- pcor.rec(x,z0,zc,method=method,na.rm=na.rm)
		ryz0.zc <- pcor.rec(y,z0,zc,method=method,na.rm=na.rm)
		
		rxy.z <- (rxy.zc - rxz0.zc*ryz0.zc)/( sqrt(1-rxz0.zc^2)*sqrt(1-ryz0.zc^2) )
		return(rxy.z)
	}			
}	
pcor <- function (y, xmat) {
	n <- dim(xmat)[2]
	resmat <- matrix(0, 1, n)
	colnames(resmat) <- colnames(xmat)
	for (i in seq(1:n)) 
		resmat[i] <- pcor.test(y, xmat[,i], xmat[,-i], method="spearman")$estimate
	attr(resmat,"N") <- dim(xmat)[1]
	return(resmat)
}
	
plot.pcor <- function (PRCC, N=attr(PRCC,"N"), nvars=dim(PRCC)[2], ...) {
	plot.pars <-list(...)
	if (is.null(plot.pars[["main"]])) {
			barplot(PRCC, ylim=c(-1,1), main="PRCC analysis", ...)
	} else {
			barplot(PRCC, ylim=c(-1,1), ...)
	}
	#Desenha as linhas a partir das quais o PRCC eh significativo
	tval <- function (prcc, df) {
		return(prcc*(sqrt(df/(1-prcc*prcc))) - qt(0.975,df) ) 
	}
	psig <- uniroot(tval, c(0,1), N-nvars-1)$root 
	abline(h=0);abline(h=psig,lty=2);abline(h=-psig,lty=2)
}
