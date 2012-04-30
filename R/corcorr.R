### Funcoes para manipulacao das variaveis visando correcao de correlacoes

s <- function (obj, x1, x2) {
	tmp <- obj[x1]
	obj[x1] <- obj[x2]
	obj[x2] <- tmp
	return (obj)
}

getE <- function (R, l, COR, i, j) {
	N <- length(R[,1])
	Tj <- matrix(0,1,length(R[1,]))
	for (m in 1:(l-1))
		Tj[,m] <- (1/N*(sum(s(R[,l],i,j)*R[,m])) - COR[l,m])^2
	return (sum(Tj))
}

# Uso: Para forcar correlacao = 0
# newvars <- LHScorcorr(vars)
# Para aproximar de uma matriz de correlacao dada:
# newvars <- LHScorcorr(vars,COR)

# NAO MEXA nos parametros l e it

LHScorcorr <- function (vars, COR = matrix(0,length(vars),length(vars)), l = 2, eps = 0.0005, it = 1) {
	# Condicao de parada: terminamos a correcao
	if (l == length(vars[1,]) + 1) {
		return (vars);
	}
	# Condicao de skip: correlacoes pequenas o suficiente:
	if (max(abs(cor(vars)[l,1:(l-1)])) < eps) {
		return (LHScorcorr (vars, COR, l = l + 1, eps = eps, it = 1));
	}
	# Condicao de parada: iteracao maxima atingida para a mesma variavel
	if (it > 20) { 
		print("ERROR: correlacao nao converge para o esperado apos numero maximo de iteracoes");
		return (LHScorcorr (vars, COR, l = l + 1, eps = eps, it = 1));
	}
	print(paste("INFO: Realizando correcao de correlacao para l =",l,"/",length(vars[1,])))
	# Aqui comeca o trabalho para corrigir as cors da var[,l]
	# Normaliza as variaveis
	N <- length(vars[,1])
	R<- matrix(0,N,length(vars[1,]))
	for (i in (1:l))
		R[,i] <- (vars[,i] - mean(vars[,i]))/sd(vars[,i])
	# Varre a matrix N por N procurando o menor erro
	# min* usados para guardar valor e posicao do minimo ateh agora
	minE <- +Inf
	mini <- 0
	minj <- 0
	for (i in (1:(N-1))) {
		for (j in ((i+1):N)) {
			E <- getE(R, l, COR, i, j)
			if (E < minE) {mini <- i; minj <- j; minE <- E}
		}
	}
	vars[,l] <- s(vars[,l], mini, minj)
	# Variavel corrigida, repete o procedimento com a mesma variavel
	# ateh a correlacao estar pequena o suficiente
	LHScorcorr(vars,COR,l=l,eps=eps,it=it+1)
}
