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

### Geradores das amostragens e funcoes acessorias:
LHSsample <- function (N, qprob=qunif, ...) { # name is required for backward compat only
	pos <- qprob((1:N)/N - 1/N/2, ...)
	pos <- sample(pos);
	attr(pos,"distribution") <- qprob # LHSextends depende desse attr
#	attr(pos,"max") <- qprob(0.99, ...)
#	attr(pos,"min") <- qprob(0.01, ...)
#	attr(pos,"name") <- name
	return (pos);
}

LHSextend <- function (old, ...) {
	qprob <- attr(old, "distribution")
	newpos <- qprob(c((1:N)/N - 5/N/6, (1:N)/N - 1/N/6), ...)
	newpos <- sample(newpos)
	attr(newpos,"distribution") <- qprob
	#     attr(newpos,"max") <- attr(old, "name")
	#     attr(newpos,"min") <- attr(old, "name")
	#     attr(newpos,"name") <- attr(old, "name")
	return (newpos);
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
oneCorPlot <- function(res, var, name) {
#		if(is.null(name)) name=attr(var,"name");
	plot(var,res, xlab=name)
	l <- lm(res~var)
	if (!is.na(coefficients(l)[2])) 
			abline((lm(res~var)))
	mtext(paste("Cor:",format(cor(var,res), digits=2), sig(res~var)))
}
corPlot <- function (vars, res, ...) {
	nl <- floor(sqrt(length(vars)))
	nc <- ceiling(length(vars)/nl)
	par(mfrow=c(nl,nc), ...)
	for (i in 1:nl) for (j in 1:nc) {
			index <- (i-1)*nc+j;
		if (index <= length(vars)) oneCorPlot(res,vars[[index]], names(vars)[index])
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
	for (i in 1:l) {
		for (j in (l+1):2) {
			if (i >= j) {
				plot.new()
			} else {
				oneTestPlot(vars[[j]],vars[[i]],test, convex, names(vars)[j], names(vars)[i])
			}
		}
	}
	par(opar)
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
	for (i in 1:l) {
		for (j in (l+1):2) {
			if (i >= j) {
				plot.new()
			} else {
				oneGradPlot(vars[[j]],vars[[i]],res, names(vars)[j], names(vars)[i])
			}
		}
	}
	par(opar)
}

# onePredPlot <- function(p1, p2, c0, c1, c2, res) {
#     p1pred <- rep(1:N/N*(limits(p1)[2]-limits(p1)[1]) + limits(p1)[1] - 1/N*(limits(p1)[2]-limits(p1)[1])/2, N)
#     p2pred <- rep(1:N/N*(limits(p2)[2]-limits(p2)[1]) + limits(p2)[1] - 1/N*(limits(p2)[2]-limits(p2)[1])/2, each=N)
#     prd <- c0 + p1pred*c1 + p2pred*c2
# # # Normalizacao
#     Max <- max(res)
#     Min <- min(res)
#     prd <- (prd-Min) / (Max-Min)
#     prd[prd < 0] <- 0
#     prd[prd > 1 ] <- 1
#     plot(0,0, xlim=limits(p1),ylim=limits(p2),xlab=paste("Valores de ",attr(p1,"name")),ylab=paste("Valores de ",attr(p2,"name"))) 
#     points(p1pred,p2pred,col=rgb(prd,0,0),pch=15)
#     points(p1,p2,col=rgb((res-Min)/(Max-Min),0,0),pch=23)
# }
# 
# predPlot <- function(vars, res, model, ...) {
#     l <- length(vars)-1
#     opar <- par(mfrow=c(l,l), ...)
#     model$coefficients -> cf
#     for (i in 1:l) {
#         for (j in (l+1):2) {
#             if (i >= j) {
#                 plot.new()
#             } else {
#                 rtmp <- res
#                 for (k in 0:l+1) { if(k != i & k != j) {
#                     rtmp <- rtmp - cf[k+1] * mean(vars[[k]])  } }
# onePredPlot(vars[[j]],vars[[i]],cf[1], cf[j+1],cf[i+1],rtmp)
#             }
#         }
#     }
#     par(opar)
# }
