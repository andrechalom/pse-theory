# EDIÇÃO FEV/2013
# Remoção de funções acessórias para facilitar o entendimento do script
# Retornar a svn 841 para funções antigas


library(sensitivity)
# Para sobol2007
genX <- function (factors, n, r, q.arg) {
	p <- length(factors)
	X <- as.data.frame(matrix(nrow = n, ncol = p))
	colnames(X) <- factors
	for (i in 1:p) {
			X[,i] <- do.call(r[i], c(n, q.arg[[i]]))
	}
	return(X)
}

mysobol <- function (model, factors, n, r, q.arg) sobol2007(model, genX(factors, n, r, q.arg), genX(factors, n, r, q.arg), nboot=1000)

myplot <- function(x, ylim = c(0, 1), ...) {
		bar.col <- c("white","orange")
		if(class(x)=="fast99") {
				if (! is.null(x$y)) {
						S <- rbind(x$D1 / x$V, 1 - x$Dt / x$V - x$D1 / x$V)
						colnames(S) <- colnames(x$X)
						barplot(S, ylim = ylim, col = bar.col)
				}
		}
		if (class(x)=="sobol2007") {
				S <- rbind(x$S$original, x$T$original - x$S$original)
				S[S<0] <- 0;
				colnames(S) <- colnames(x$X)
				b<-barplot(S, ylim= ylim, col=bar.col)
				smin <- x$S$"min. c.i."
				smax <- x$S$"max. c.i."
				tmin <- x$T$"min. c.i."
				tmax <- x$T$"max. c.i."
				for (i in 1:length(colnames(S))) {
						lower <- min(smin[i], tmin[i])
						upper <- max(smax[i], tmax[i])
						arrows(b[i],lower, b[i], upper, angle=90, code=3, length=0.1) #, lty=2)
						#				lower <- min(max(smin[i], tmin[i]), x$S$original[i])
						#				upper <- max(min(smax[i], tmax[i]), x$S$original[i])
						#				segments(b[i],lower, b[i], upper, lwd=2)
				}
		}
		legend("topright", c("main effect", "interactions"), fill = bar.col)
}

# Para poder plotar PRCC com mais liberdade (cor, titulos, etc)

nodeplot <- function(x, xlim = NULL, ylim = NULL, labels = TRUE,
					 col = par("col"), pch = 21, bg = "white",
					 add = FALSE, at = NULL, ...) {
		n <- nrow(x)
		if (is.null(xlim)) {
				xlim <- c(1, n)
		}
		if (is.null(ylim)) {
				ylim <- c(min(x), max(x))
		}
		if (is.null(at)) {
				at <- 1 : n
		}
		if (add) {
				par(new = TRUE)
		}

		# axes

		plot(0, xlim = xlim, ylim = ylim, axes = FALSE,
			 xlab = "", ylab = "", type = "n", ...)
		if (class(labels) == "logical") {
				if (labels) {
						axis(side = 1, at = at, labels = rownames(x))
				} else {
						axis(side = 1, at = at, labels = FALSE, tick = FALSE)
				}
		} else if (class(labels) == "character") {
				axis(side = 1, at = at, labels = labels)
		}
		axis(side = 2)
		box()

		# bias

		if ("bias" %in% colnames(x)) {
				xx <- x[["original"]] - x[["bias"]]
		} else {
				xx <- x[["original"]]
		}

		# confidence intervals

		if (("min. c.i." %in% colnames(x)) & "max. c.i." %in% colnames(x)) {
				for (i in 1 : n) {
						lines(c(at[i], at[i]), c(x[["min. c.i."]][i], x[["max. c.i."]][i]),
							  col = col)
				}
		}

		# points

		points(at, xx, col = col, pch = pch, bg = bg)
}

########################CORCORR

# Uso: Para forcar correlacao = 0
# newvars <- LHScorcorr(vars)
# vars deve ser uma data.frame ou objeto que pode ser coerced para data.frame
# O valor de retorno eh sempre um data.frame com os mesmos rownames de vars
# Para aproximar de uma matriz de correlacao dada:
# newvars <- LHScorcorr(vars,COR)

system("R CMD SHLIB corcorr.c")

# NAO MEXA nos parametros l e it
LHScorcorr <- function (vars, COR = matrix(0,dim(vars)[2],dim(vars)[2]), l = 2, eps = 0.0005, it = 1, echo=FALSE, maxIt = 2*sqrt(dim(vars)[1])) {
		N <- dim(vars)[1]; M <- dim(vars)[2]
		my.names <- names(vars)
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
		vars <- as.data.frame(matrix(V$vars, nrow=N, ncol=M))
		names(vars) <- my.names
		if (V$FLAGSTOP == 1) { # O procedimento convergiu
				return(LHScorcorr (vars, COR, l = l + 1, eps = eps, it = 1, echo=echo));
		} else {
				# Variavel corrigida, repete o procedimento com a mesma variavel
				# ateh a correlacao estar pequena o suficiente
				LHScorcorr(vars,COR,l=l,eps=eps,it=it+1, echo=echo)
		}
}
###################################

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
											  N="numeric", data="data.frame", factors="character", q="character", 
											  q.arg="list", COR="matrix", eps="numeric", model="function", res="numeric",
											  prcc = "list")
)

results <- function(LHS) { return (LHS@res) }
setMethod(
		  f="pcc",
		  signature("LHS"),
		  definition = function (X, y=NULL, rank=TRUE,nboot=1000, conf=0.95) {
				  require(sensitivity)
				  iprcc <- sensitivity::pcc(X@data, as.vector(X@res), nboot=nboot, rank=rank, conf=conf)
				  return(iprcc);
		  }
		  )

setMethod(
		  f="ecdf",
		  signature("LHS"),
		  definition = function (x) {
				  return(ecdf(x@res));
		  }
		  )

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

LHS <- function (model, factors, N, q, q.arg, COR=matrix(0, length(factors), length(factors)), eps=0.0005, nboot=0) {
		L <- as.data.frame(matrix(nrow=N, ncol=length(factors)));
		names(L) <- factors
		for (i in 1:length(factors)) 
				L[,i] <- sample(do.call(q[[i]], c(list(p = 1:N/N-1/N/2), q.arg[[i]])))
		L <- LHScorcorr(L, COR,eps=eps); 
		res <- model(L);
		prcc <- pcc(L, as.vector(res), rank=TRUE);
		class(prcc) <- "list";
		X <- new(Class="LHS", N=N, data=L,factors=factors, q=q, q.arg=q.arg, COR=COR, eps=eps, model=model, res=res, prcc=prcc);
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

