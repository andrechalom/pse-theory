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

# Dada uma amostragem gerada por sample*, devolve os limites de plotagem para a variavel
limits <- function (pos) {
	if (attr(pos,"distribution") == "uniform") {
		return (c(attr(pos,"min"),attr(pos,"max")))
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
	par(mfrow=c(l,l), ...)
	for (i in 1:l) {
		for (j in (l+1):2) {
			if (i >= j) {
				plot.new()
			} else {
				oneTestPlot(vars[[j]],vars[[i]],test)
			}
		}
	}
}
