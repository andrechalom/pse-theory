# EDIÇÃO FEV/2013
# Remoção de funções acessórias para facilitar o entendimento do script
# Retornar a svn 841 para funções antigas

# Edicao MAR/2013 - funcionalidades para que a funcao "model" possa retornar
# um array
library(sensitivity)
source("corcorr.R")
source("sobol.R")
# Dada uma amostragem gerada por sample*, devolve os limites de plotagem para a variavel
limits <- function (pos) {
		return (c(min(pos), max(pos)))
}

### Plots e analises
sig <- function (form) {
		l <- lm(form)
		if (is.na(coefficients(l)[2])) return (" ");
		f <- anova(lm(form))$"Pr(>F)"[1]
		if (f < 0.001) return ("***");
		if (f < 0.01) return ("**");
		if (f < 0.05) return ("*");
		if (f < 0.1) return (".");
		return (" ");
}
oneCorPlot <- function(res, var, name, log, first, ...) {
		# TODO: Melhorar essa funcao que deforma o primeiro da linha!!
		if (first) par(mar=c(5,2,4,0.5), yaxt='s')
		else par (mar=c(5,0.5,4,0.5), yaxt='n')
		plot(var,res, xlab=name, log=log, ylab= "", ...)
		l <- lm(res~var)
		if (!is.na(coefficients(l)[2])) 
				abline((lm(res~var)))
		mtext(paste("Cor:",format(cor(var,res), digits=2), sig(res~var)))
}
corPlot <- function (vars, res=NULL, log="", ...) {

		if (class(vars)=="LHS") return(corPlot(vars@data, vars@res, log=log, ...))

		nl <- floor(sqrt(length(vars)))
		nc <- ceiling(length(vars)/nl)
		opar <- par(mfrow=c(nl,nc), pch='.', ...)
		on.exit(par(opar))
		for (i in 1:nl) for (j in 1:nc) {
				index <- (i-1)*nc+j;
				if (index <= length(vars)) oneCorPlot(res,vars[[index]], names(vars)[index], log, first = (j == 1), ...)
		}
}

# Funcao acessoria - nao eh para uso externo
oneTestPlot <- function(p1, p2, test, convex=FALSE, p1name, p2name) {
		if (convex) {
				require(spatstat)
				W <- owin(c(min(p1[test]),max(p1[test])), c(min(p2[test]),max(p2[test])))
				P <- ppp(p1[test],p2[test],window=W)
				cTrue <-convexhull(P)
				W <- owin(c(min(p1[!test]),max(p1[!test])), c(min(p2[!test]),max(p2[!test])))
				P <- ppp(p1[!test],p2[!test],window=W)
				cFalse <- convexhull(P)
		}
		plot(0,0,
			 xlim=limits(p1),ylim=limits(p2),
			 xlab=paste("Valores de ",p1name),ylab=paste("Valores de ",p2name), type='n')
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
		on.exit(par(opar))
		for (i in 1:l) {
				for (j in (l+1):2) {
						if (i >= j) {
								plot.new()
						} else {
								oneTestPlot(vars[[j]],vars[[i]],test, convex, names(vars)[j], names(vars)[i])
						}
				}
		}
}

oneGradPlot <- function (p1, p2, res, p1name, p2name) {
		Max <- max(res)
		Min <- min(res)
		plot(0,0,
			 xlim=limits(p1),ylim=limits(p2),
			 xlab=paste("Valores de ",p1name),ylab=paste("Valores de ",p2name))
		points(p1,p2,col=rgb((res-Min)/(Max-Min),0,0),pch=19)
}

# Plota o resultado da simulacao pelo hipercubo
gradPlot <- function(vars, res, ...) {
		l <- length(vars)-1
		opar <- par(mfrow=c(l,l), ...)
		on.exit(par(opar))
		for (i in 1:l) {
				for (j in (l+1):2) {
						if (i >= j) {
								plot.new()
						} else {
								oneGradPlot(vars[[j]],vars[[i]],res, names(vars)[j], names(vars)[i])
						}
				}
		}
}


sbma <- function (sample1, sample2, absolute=TRUE) {
		if(class(sample1) == "fast99") 
				return(data.frame(Di=sbma(sample1$D1, sample2$D1), Dt=sbma(sample1$Dt, sample2$Dt)))
		if(class(sample1) == "sobol2007") 
				return(data.frame(S=sbma(sample1$S$original, sample2$S$original), T=sbma(sample1$T$original, sample2$T$original)))
		if(class(sample1) == "LHS") {
				x1 <- sample1@prcc$PRCC$original; 
				x2 <- sample2@prcc$PRCC$original;
		}
		else {x1 <- sample1; x2<-sample2;}
		if (absolute) {x1 <- abs(x1); x2 <- abs(x2);}
		if (length(x1) != length(x2)) 
				stop("Sample sizes must be the same!");
		n <- length(x1);
		R <- rank(x1);
		S <- rank(x2);
		v1 <- -(4*n+5)/(n-1);
		v2 <- 6/(n^3-n);
		sbc <- v1+v2*(sum(R*S*(4-(R+S)/(n+1))));
		return (sbc)
}

setClass("LHS", representation=representation(
											  N="numeric", data="matrix", factors="character", q="character", 
											  q.arg="list", COR="matrix", eps="numeric", model="function", res="matrix",
											  res.names="character", prcc = "list")
)

get.results <- function(LHS) { return (LHS@res) }
get.data <- function(LHS) {return (LHS@data) }
setMethod(
		  f="pcc",
		  signature("LHS"),
		  definition = function (X, y=NULL, rank=TRUE,nboot=1000, conf=0.95) {
				  require(sensitivity)
				  iprcc <- sensitivity::pcc(X@data, as.vector(X@res), nboot=nboot, rank=rank, conf=conf)
				  return(iprcc);
		  }
		  )

# setMethod(
#           f="ecdf",
#           signature("LHS"),
#           definition = function (x) {
#                   return(ecdf(x@res));
#           }
#           )


plotecdf <- function (LHS, ...) {
		nres <- dim(LHS@res)[2]
		nl <- floor(sqrt(nres))
		nc <- ceiling(nres/nl)
		opar <- par(mfrow=c(nl,nc), pch='.', ...)
				min <- min(LHS@res)
				max <- max(LHS@res)
	for(i in 1:nres) 
			plot(ecdf(LHS@res[,i]), xlim=c(min, max), main='', xlab=LHS@res.names[i])
}

setMethod(
		  f="kruskal.test",
		  signature("LHS"),
		  definition = function (x, Ncats=10, do.plot=TRUE, ...) {
				  k <- NA;
				  df <- Ncats-1;
				  for (i in 1:length(x@data)) {
						  cats <- cut(x@data[,i], breaks=quantile(x@data[,i], seq(0,1,1/Ncats)))
						  k[i]<-(kruskal.test(x@res~cats))$statistic;
				  }
				  names(k) <- x@factors;
				  if(do.plot) {
						  plot(k, ylim=c(0, max(k)*1.1), xaxt='n')
						  axis(1, at=1:length(x@data), labels=x@factors)
						  abline(h=qchisq(0.05, df=df), lty=2)
				  }
				  return (k);
		  }
		  )

LHS <- function (model, factors, N, q, q.arg, res.names, COR=matrix(0, length(factors), length(factors)), eps=0.0005, nboot=0) {
		L <- matrix(nrow=N, ncol=length(factors));
		colnames(L) <- factors
		for (i in 1:length(factors)) 
				L[,i] <- sample(do.call(q[[i]], c(list(p = 1:N/N-1/N/2), q.arg[[i]])))
		L <- LHScorcorr(L, COR,eps=eps); 
		res <- t(model(L));
		#prcc <- pcc(L, res, rank=TRUE);
		prcc <- "to be implemented"
		class(prcc) <- "list";
		X <- new(Class="LHS", N=N, data=L,factors=factors, q=q, q.arg=q.arg, COR=COR, eps=eps, model=model, res=res, prcc=prcc, res.names=res.names);
		return(X);
}

LHSextend<- function (oldLHS) {
		newL <- as.data.frame(matrix(nrow=2*oldLHS@N, ncol=length(oldLHS@factors)));
		names(newL) <- oldLHS@factors
		for (i in 1:length(oldLHS@factors)) 
				newL[,i] <- do.call(oldLHS@q[[i]], c(list(p =c((1:N)/N - 5/N/6, (1:N)/N - 1/N/6)), oldLHS@q.arg[[i]]))
		newL <- LHScorcorr(newL, COR=oldLHS@COR, eps=oldLHS@eps);
		newres <- oldLHS@model(newL);
		oldLHS@N <- oldLHS@N*3;
		oldLHS@data <- rbind(oldLHS@data, newL);
		oldLHS@res <- c(oldLHS@res, newres);
		prcc <- pcc(oldLHS@data, as.vector(oldLHS@res), rank=TRUE)
		class(prcc) <- "list";
		oldLHS@prcc <- prcc
		return(oldLHS);
}

