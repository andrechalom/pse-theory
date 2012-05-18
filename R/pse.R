### Geradores das amostragens e funcoes acessorias:
LHSsample <- function (N, name, qprob, ...) {
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

# Funcao de correcao de correlacoes
source("corcorr.R")

# Metodo de calcular a PRCC:
source("pcor.R")
pcor <- function (y, xmat) {
	n <- dim(xmat)[2]
	resmat <- matrix(0, 1, n)
	colnames(resmat) <- colnames(xmat)
	for (i in seq(1:n)) 
		resmat[i] <- pcor.test(y, xmat[,i], xmat[,-i], method="spearman")$estimate
	return(resmat)
}
	
plot.pcor <- function (PRCC, N, nvars) {
	barplot(PRCC, ylim=c(-1,1), main="PRCC analysis")
	#Desenha as linhas a partir das quais o PRCC eh significativo
	tval <- function (prcc, df) {
		return(prcc*(sqrt(df/(1-prcc*prcc))) - qt(0.975,df) ) 
	}
	psig <- uniroot(tval, c(0,1), N-nvars-1)$root ### df deve ser alterado NA MAO
	abline(h=0);abline(h=psig,lty=2);abline(h=-psig,lty=2)
}
