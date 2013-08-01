# EDIÇÃO FEV/2013
# Remoção de funções acessórias para facilitar o entendimento do script
# Retornar a svn 841 para funções antigas

# Edicao MAR/2013 - funcionalidades para que a funcao "model" possa retornar
# um array
library(sensitivity)
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
		my.names <- colnames(vars)
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
### Plots e analises
oneCorPlot <- function(res, var, name, log, add.lm, ...) {
		par (mar=c(4,4,2,0.5))
		plot(var,res, xlab=name[1], log=log, ylab= name[2], ...)
		if (add.lm) {
			l <- lm(res~var)
			if (!is.na(coefficients(l)[2])) 
					abline((lm(res~var)))
		}
}
corPlot <- function (vars, res=NULL, index.data=NULL, index.res=NULL,log="", add.lm=TRUE, ...) {
		if (class(vars)=="LHS") {
				if (is.null(index.data)) index.data <- 1:dim(vars@data)[2]
				if (is.null(index.res)) index.res <- 1:dim(vars@res)[2]
				return(corPlot(vars@data[index.data], vars@res[index.res], log=log, add.lm=add.lm, ...))
		}
		res <- as.data.frame(res)
		if (is.null(dim(vars))) {nplots <-1} else {nplots <- dim(vars)[2]}
		if (is.null(dim(res))) {nplots <-nplots*1} else {nplots <- nplots* dim(res)[2]}

		nl <- floor(sqrt(nplots))
		nc <- ceiling(nplots/nl)
		opar <- par(mfrow=c(nl,nc), pch='.', ...)
		on.exit(par(opar))

		index.var <- 1
		index.res <- 1
		for (i in 1:nl) for (j in 1:nc) {
				oneCorPlot(res[,index.res],vars[, index.var], c(names(vars)[index.var],names(res)[index.res]), log, add.lm, ...)
				index.res <- index.res +1
				if (index.res > dim(res)[2]) {index.res <- 1; index.var <- index.var + 1}
				if (index.var > dim(vars)[2]) return();
		}
}

#                       Nodeplot: anti-boxplot
#                         Gilles Pujol 2006
#   Editado Chalom 2013 para bugfix na versao 1.6-1, aguardando updates

plot.pcc <- function(x, ylim = c(-1,1), main=(if ("PCC" %in% names(x)) "PCC" else "PRCC"), ...) {  
		if ("PCC" %in% names(x)) {
				nodeplot(x$PCC, ylim = ylim, main=main, ...)
		}else if ("PRCC" %in% names(x)) {
				nodeplot(x$PRCC, ylim = ylim, main=main, ...)
		}  
}

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
       xlab = "", type = "n", ...)
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

sbma <- function (sample1, sample2, absolute=TRUE) {
		if(class(sample1) == "LHS") {
				sb <- array();
				for (i in 1:dim(sample1@res)[2])
						sb[i] <- sbma(sample1@prcc[[i]]$PRCC$original, sample2@prcc[[i]]$PRCC$original);
				return (sb);
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
											  q.arg="list", COR="matrix", eps="numeric", model="function", res="data.frame",
											  res.names="character", prcc = "list")
)

get.results <- function(LHS) { return (LHS@res) }
get.data <- function(LHS) {return (LHS@data) }

setMethod(
		  f="pcc",
		  signature("LHS"),
		  definition = function (X, y=NULL, rank=TRUE,nboot=1000, conf=0.95) {
				  require(sensitivity)
				  f <- function(res) sensitivity::pcc(X@data, as.vector(res), nboot=nboot, rank=rank, conf=conf)
				  return(apply(X@res, 2, f))
		  }
		  )

plotecdf <- function (LHS, stack=FALSE, index.res =1:dim(LHS@res)[2], col=index.res, ...) {
		require(Hmisc)
		require(ks)
		if (stack) {
			dat <- vec(LHS@res[,index.res])
			g <- rep(index.res, each=dim(LHS@res)[1])
			Ecdf(dat, group=g, col=col)
		} else Ecdf(LHS@res[,index.res])
}

plotprcc <- function (LHS, bg ='orange', index.res = 1:dim(LHS@res)[2], ...) {
		nres <- length(index.res)
		nl <- floor(sqrt(nres))
		nc <- ceiling(nres/nl)
		opar <- par(mfrow=c(nl,nc), ...)
		on.exit(par(opar))
	for(i in 1:nres) 
	{
			plot.pcc(LHS@prcc[[index.res[i]]], bg=bg, main=LHS@res.names[[index.res[i]]])
			abline(h=0, lty=2)
	}
}

LHS <- function (model, factors, N, q, q.arg, res.names=NULL, COR=matrix(0, length(factors), length(factors)), eps=0.0005, nboot=0) {
		L <- as.data.frame(matrix(nrow=N, ncol=length(factors)));
		colnames(L) <- factors
		for (i in 1:length(factors)) 
				L[,i] <- sample(do.call(q[[i]], c(list(p = 1:N/N-1/N/2), q.arg[[i]])))
		L <- LHScorcorr(L, COR,eps=eps); 
		res <- t(model(L));
		if (dim(res)[1] == 1) res <- t(res); # Caso seja retornado um unico valor
		res <- as.data.frame(res);
		if(!is.null(res.names))
			colnames(res) <- res.names
		else res.names <- paste("V", 1:dim(res)[2], sep="")
		# Calculates the PRCC for each response variable
		f <- function(r) sensitivity::pcc(L, r, nboot=nboot, rank=T)
		prcc <- apply(res, 2, f);
		###
		X <- new(Class="LHS", N=N, data=L,factors=factors, q=q, q.arg=q.arg, COR=COR, eps=eps, model=model, res=res, prcc=prcc, res.names=res.names);
		return(X);
}

target.sbma <- function(target, model, factors,  q, q.arg, res.names=NULL, COR=matrix(0, length(factors), length(factors)), eps=0.0005,init=length(factors)+2, inc=100, FUN=min) {
		#initial LHS
		N = init
		print("INFO: initial run...")
		oldL <- LHS(model, factors, N, q, q.arg, res.names, COR, eps, nboot=0)
		while (TRUE) {
				N = N + inc
				print(paste("INFO: LHS with N =", N));
				newL <- LHS(model, factors, N, q, q.arg, res.names, COR, eps, nboot=0)
				s <- FUN(sbma(newL, oldL))
				print(paste("sbma of ", round(s,3)," (target ",target,")", sep=""))
				if (s >= target) return (newL);
				oldL <- newL;
		}
}
						
# Sensitivity and elasticity (Caswell)
estim.slr <- function(data, i = 1:nrow(data)) {  
  d <- data[i, ]
  p <- ncol(d) - 1
  slr <- numeric(p)
  for (j in 1:p) {
    Xtildej.lab <- paste(colnames(d)[c(-1, -j-1)], collapse = "+")
    lm.Y <- lm(formula(paste(colnames(d)[1], "~", Xtildej.lab)), data = d)
    lm.Xj <- lm(formula(paste(colnames(d)[j+1], "~", Xtildej.lab)), data = d)
	y = d[1] - fitted(lm.Y)
	x = d[j+1] - fitted(lm.Xj)
	slr[j] <- coef(lm(y[,1] ~ x[,1]))[2]
  }
  return(slr)
}

slr <- function (LHS, nboot = 0, conf=0.95) {
	data <- cbind(Y=get.results(LHS), get.data(LHS))
	if (nboot == 0) {
		slr <- data.frame(original = estim.slr(data))
		rownames(slr) <- colnames(get.data(LHS))
	} else {
		boot.slr <- boot(data, estim.slr, R = nboot)
		slr <- bootstats(boot.slr, conf, "basic")
		rownames(slr) <- colnames(get.data(LHS))
	}
	out <- list(X = get.data(LHS), y = get.results(LHS), nboot = nboot, conf = conf,
				call = match.call(), slr = slr)
	class(out) <- "slr"
	return(out)
}

elast <- function (LHS) {
	s <- slr(LHS)
	f <- apply (get.data(LHS), 2, mean, na.rm=T)
	m <- apply(get.results(LHS),2, mean, na.rm=T)
	return(s$slr$original * f / m)
}

plotelast <- function (LHS) {
	s <- slr(LHS)
	f <- apply (get.data(LHS), 2, mean, na.rm=T)
	m <- apply(get.results(LHS),2, mean, na.rm=T)
	barplot(s$slr$original * f / m, ylim=c(-1,1))
}
