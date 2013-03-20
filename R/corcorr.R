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
		vars <- matrix(V$vars, nrow=N, ncol=M)
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
