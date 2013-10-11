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
	if (nboot == 0) {
		slr <- list()
		for(i in 1:dim(get.results(LHS))[2]) {
			data <- cbind(Y = get.results(LHS)[i], get.data(LHS))
			slr[[i]] <- estim.slr(data)
		}
		slr <- as.data.frame(slr)
		rownames(slr) <- colnames(get.data(LHS))
		colnames(slr) <- colnames(get.results(LHS))
	} else {
		data <- cbind(Y=get.results(LHS), get.data(LHS))
		## NAO funciona para mais de uma variavel resposta
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
	barplot(s$slr * f / m, ylim=c(-1,1))
}
