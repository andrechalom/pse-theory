# EDIÇÃO FEV/2013
# Remoção de funções acessórias para facilitar o entendimento do script
# Retornar a svn 841 para funções antigas

# Edicao MAR/2013 - funcionalidades para que a funcao "model" possa retornar
# um array

# Edicao out/2013 - pse agora eh um pacote, restam aqui funcoes tentativas
# e não testadas
library(pse)
########################CORCORR

# Para sobol2007
library(sensitivity)
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

# Para PRCC

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


# Sensitivity and elasticity (Caswell)
### WORK IN PROGRESS: adaptar as funcoes para multiplas respostas
elast <- function (LHS) {
	s <- pic(LHS)
	f <- apply (get.data(LHS), 2, mean, na.rm=T)
	m <- apply(get.results(LHS),2, mean, na.rm=T)
	return(s[[1]]$pic$original * f / m)
}

plotelast <- function (LHS) {
	s <- pic(LHS)
	f <- apply (get.data(LHS), 2, mean, na.rm=T)
	m <- apply(get.results(LHS),2, mean, na.rm=T)
	barplot((s[[1]]$pic * f / m)[,1], ylim=c(-1,1))
}


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

