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

