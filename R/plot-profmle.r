#=======================================================================
plot.profmle <- function(mleobj, nseg=20, ratio=log(8), which=NULL, auto.mfrow=TRUE, ... )
{
	# Guarda os nossos parametros graficos atuais e garante que eles
	# irao retornar apos a saida da funcao:
	opar <- par(no.readonly=TRUE)
	on.exit(par(opar))

	if( class(mleobj)[1] != "profile.mle" &
	   class(mleobj)[1] != "profile.mle2") 
			stop( "Object should have class \'profile.mle\' or \'profile.mle2\'")
	mleprof <- mleobj@profile
	npar <- length(mleprof)
	if( is.null(which) )
	{
			parseq = 1:npar
			if (auto.mfrow) {
				nl <- floor(sqrt(npar))
				nc <- ceiling(npar/nl)
				par(mfrow=c(nl,nc))
			}
	}
	else
			parseq = which
	for(i in parseq)
	{
			tmp <- mleprof[i][[1]]
		varname <- names(mleprof[i])
			y <- tmp[,1]^2/2
			x <- (tmp[,2][,i])
			interpol = spline(x, y, n=nseg*length(x) )
			plot(interpol, 
				 type="l", 
				 xlab=varname, 
				 ylab="Log-Verossimilhança Negativa Relativa",
				 col="red",
				 ...
				 )

			# Codigo antigo para tracar o intervalo de verossimilhanca
			#         x.verint = interpol$x[ interpol$y <= ratio ]
			#         y.verint = rep( ratio, length(x.verint) ) 
			#         lines(x.verint, y.verint , col="blue", lty=2)
			#         lines(rep(x.verint[1],2), c(-1,ratio), col="blue", lty=2)
			#         lines(rep(x.verint[length(x.verint)],2), c(-1,ratio), col="blue", lty=2)

			# Determinamos em que intervalos nossa linha de y = ratio 
			# intercepta o spline
			l <- length(interpol$y)
			change <- (interpol$y - ratio)[2:l] * (interpol$y - ratio)[1:(l-1)]
			endpoints <- which(change < 0)
			# interpol$x[change] guarda o limite INFERIOR de cada intervalo,
			# mas vamos tracar as linhas no ponto medio de cada intervalo:
			# Talvez haja maneiras mais elegantes de fazer isso??
			corr <- (interpol$x[2]-interpol$x[1])/2
			# Para cada intervalo encontrado, desenhamos os limites do intervalo de verossimilhanca
			if (length(endpoints)) {
				for (j in 1:(length(endpoints)/2)) {
						lower <-interpol$x[endpoints[(2*j)-1]]+corr
						upper <- interpol$x[endpoints[2*j]]+corr
						lines(c(lower,upper ),c(ratio, ratio), col="blue", lty=2)
						lines(rep(lower,2), c(-1, ratio), col="blue", lty=2)
						lines(rep(upper,2), c(-1, ratio), col="blue", lty=2)
				}
			}
	}
}
