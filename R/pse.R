# EDIÇÃO FEV/2013
# Remoção de funções acessórias para facilitar o entendimento do script
# Retornar a svn 841 para funções antigas

# Edicao MAR/2013 - funcionalidades para que a funcao "model" possa retornar
# um array

# Edicao out/2013 - pse agora eh um pacote, restam aqui funcoes tentativas
# e não testadas
library(pse)
########################CORCORR

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

