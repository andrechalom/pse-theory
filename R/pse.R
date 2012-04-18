### Geradores das amostragens e funcoes acessorias:

# sampleunif: gera uma amostra proveniente de N intervalos de uma
# distribuicao uniforme
# Devolve os valores como array, e seta os atributos distribution, max e min
sampleunif <- function (N, xmin, xmax, name) {
	limits <- 0:N/N*(xmax-xmin) + xmin
	pos <- array(0)
	for (i in 1:N) pos[i] <- runif(1,limits[i],limits[i+1])
	pos <- sample(pos);
	attr(pos,"distribution") <- "uniform"
	attr(pos,"max") <- xmax
	attr(pos,"min") <- xmin
	attr(pos,"name") <- name
	return (pos);
}

samplenorm <- function (N, mean, sd, name) {
	pos <- qnorm((1:N)/N - 1/N/2, mean, sd )
	pos <- sample(pos);
	attr(pos,"distribution") <- "normal"
	attr(pos,"mean") <- mean
	attr(pos,"sd") <- sd
	attr(pos,"name") <- name
	return (pos);
}

# Dada uma amostragem gerada por sample*, devolve os limites de plotagem para a variavel
limits <- function (pos) {
	if (attr(pos,"distribution") == "uniform") {
		return (c(attr(pos,"min"),attr(pos,"max")))
	} else if (attr(pos,"distribution") == "normal") {
		m<-attr(pos,"mean"); s <- attr(pos,"sd");
		return (c(qnorm(0.01,m,s),qnorm(0.99,m,s)))
	}
	return (NULL);
}

### Plots e analises

# Funcao acessoria - nao eh para uso externo
oneTestPlot <- function(p1, p2, test) {
	plot(0,0, xlim=limits(p1),ylim=limits(p2),xlab=paste("Valores de ",attr(p1,"name")),ylab=paste("Valores de ",attr(p2,"name")))
	points(p1[test],p2[test], pch='+')
	points(p1[!test],p2[!test], pch='-')
}
# Plota o resultado de um teste logico pelo hipercubo
testPlot <- function(vars, test, ...) {
	l <- length(vars)-1
	opar <- par(mfrow=c(l,l), ...)
	for (i in 1:l) {
		for (j in (l+1):2) {
			if (i >= j) {
				plot.new()
			} else {
				oneTestPlot(vars[[j]],vars[[i]],test)
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
			#	for (k in 0:l+1) { if(k != i & k != j) {
			#	       print (paste("Corrigindo plot de ",attr(vars[[j]],"name")," x ", attr(vars[[i]],"name")," pela media de ",attr(vars[[k]],"name")))	
			#		rtmp <- rtmp - cf[k] * vars[[k]]  } }
onePredPlot(vars[[j]],vars[[i]],cf[1], cf[j+1],cf[i+1],rtmp)
			}
		}
	}
	par(opar)
}
