source("pse.R")
modelRun <- function (Xo, r, K, Time) {
	X <- Xo
	for (i in 0:Time) {
		X <- r*X[length(X)]*(1-X[length(X)]/K)
	}
	return (X)
}
# First sampling:
N <-20
quant <- (1:N)/N - 1/N/2
quant
# Triple N
newquant <- c((1:N)/N - 5/N/6, (1:N)/N - 1/N/6  )
# Construo o cubo pequeno
r <- qnorm(quant, 0, 1); r<-sample(r)
k <- qt(quant, 4, 1); k<-sample(k)
cor(r,k)
newvars <- LHScorcorr (cbind(r,k))
newvars[,1] -> r; newvars[,2] ->k;
cor(r, k)
# EXTENDENDO o cubo
newr <-qnorm(newquant, 0, 1); newr<-sample(newr)
newk <- qt(newquant, 4,1); newk<-sample(newk)
newvars <-  LHScorcorr (cbind(newr,newk))
newr <- newvars[,1]; newk<- newvars[,2];
print(cor(r,k)^2); print(cor(newr, newk)^2); print(cor(c(r,newr),c(k,newk))^2)
print(cor(c(r,newr),c(k,newk))^2 < (cor(r,k)^2 + cor(newr, newk)^2) )
mean(r); mean(newr); mean(c(newr,r))
sd(r); sd(newr); sd(c(newr,r))





### Descrevendo o espaco de parametros e gerando o hipercubo
# Definicao do numero total de amostras:
N <- 150
# r eh uniformemente distribuido entre 0 e 2
# os outros parametros seguem logica semelhante
r <- LHSsample(N, "r", qunif, 0.25, 2)
k <- LHSsample(N, "k", qunif, 10, 50)
x <- LHSsample(N, "x", qunif, 1, 10)
# O tempo final eh distribuido com forma normal de media = 100 e sd = 20
Time <- LHSsample(N, "T", qnorm, 100, 40)

apply((cbind(r,k,Time,x)), 2, mean)
apply((cbind(r,k,Time,x)), 2, sd)

# Correlacao entre as variaveis geradas. Veja que ela nao eh nec. 0
cor(cbind(r,k,Time,x))
# Correcao das variaveis para apresentarem correlacao 0
 newvars <- LHScorcorr (cbind(r,k,Time,x))
# Correcao das variaveis para r e k terem correlacao negativa
MyCor <- matrix(c(   1,-0.5,   0,   0,
		  -0.5,   1,   0,   0,
		     0,   0,   1,   0,
		     0,   0,   0,   1),4,4)
#newvars <- LHScorcorr (cbind(r,k,Time,x), MyCor)
# Uso a notacao [1:N] para manter os atributos - **TODO** refazer pra ficar menos pentelho
k[1:N] <- newvars[,2]
Time[1:N] <- newvars[,3]
x[1:N] <- newvars[,4]
print(cor(cbind(r,k,Time,x)))

### Rodando o modelo para os diferentes parametros
# O objeto res contem a saida do modelo, no nosso caso, a populacao final
res <- mapply(modelRun, x, r, k, Time)

#save(res,x,r,k,Time, file="LHS4500.Rdata")
### Plots e analises
# Analise da correlacao entre o resultado e cada variavel:
# metodos de Pearson e Spearman e porcentagem da variancia de res
# explicada por efeitos lineares e monotonicos dos inputs
cor(cbind(res,r,k,Time,x))
cor(res,r)^2 + cor(res,k)^2 + cor(res,Time)^2 + cor(res,x)^2
cor(cbind(res,r,k,Time,x), method="spearman")
cor(res,r, method="spearman")^2 + cor(res,k, method="spearman")^2 + cor(res,Time, method="spearman")^2 + cor(res,x, method="spearman")^2


# Distribuicao de probabilidades de y = f(x)
hist(res[res>0.001])
plot(density(res[res>0.0001]))
plot(ecdf(res[res>0.0001]), do.points=FALSE)
# TODO: Gerar uma estrutura mais generica para plotar graficos
# Plot da contribuicao de cada variavel, independentemente
# Os indicadores de significancia nao podem ser levados a serio demais...
corPlot(list(r,k,Time,x),res)

# Plot de sobrevivencia
# Nesta funcao, testamos uma condicao da saida do modelo: populacao
# final maior do que 0 (mais um errinho numerico)
# Note que as funcoes de plot aceitam parametros graficos com algumas
# condicoes (mais sobre isso no manual)
# A funcao plota + onde o teste eh verdadeiro, - onde eh falso
testPlot(list(r,x,Time,k), test=res>0.0001, mar=c(4,4,1,2))

# Plot da populacao final
# TODO: deixar maior controle sobre a cor que deve ser plotada, etc
gradPlot(list(r,x,Time,k), res, mar=c(4,4,1,2))

# Geramos um modelo linear para prever a populacao em qualquer ponto
# O resultado comum para N suficiente eh um modelo com efeito
# significativo e importante de r e k, e efeito desprezivel de T e x
modellm <- lm (res ~ r + x + Time + k)
print (modellm)
# O sumario do modelo traz as estatisticas da ANOVA correspondente
print (summary(modellm))
# anova(modellm)
# Plot da populacao prevista pelo modelo e realizada nas simulacoes
predPlot(list(r,x,Time,k), res, modellm)

# Mais um modelo linear, agora restrito a regiao r > 1
modellm <- lm (res ~ r + x + Time + k, subset=(r>1))
print (modellm)
print (summary(modellm))
predPlot(list(r,x,Time,k), res, modellm)

# PCC: Partial correlation coefficient between res & r
cor(
    residuals(hatr <- lm(r~k+Time+x)),
    residuals(hatres <- lm(res~k+Time+x))
)

# Rank transformations:
resr <- rank(res); rr <- rank(r); kr <- rank(k); Timer <- rank(Time); xr <- rank(x)
PRCC<-matrix(0,1,4)
colnames(PRCC)<-c("r","k","Time","x")
#PRCC between res & r
PRCC[1]<-
cor(
    residuals(lm(rr~kr+Timer+xr)),
    residuals(lm(resr~kr+Timer+xr))
)
#PRCC between res & k
PRCC[2] <-
cor(
    residuals(lm(kr~rr+Timer+xr)),
    residuals(lm(resr~rr+Timer+xr))
)
#PRCC between res & Time
PRCC[3] <-
cor(
    residuals(lm(Timer~kr+rr+xr)),
    residuals(lm(resr~kr+rr+xr))
)
#PRCC between res & x
PRCC[4] <-
cor(
    residuals(lm(xr~kr+Timer+rr)),
    residuals(lm(resr~kr+Timer+rr))
)
print(PRCC)
barplot(PRCC, ylim=c(-1,1), main="PRCC analysis")
#Desenha as linhas a partir das quais o PRCC eh significativo
tval <- function (prcc, df) {
	return(prcc*(sqrt(df/(1-prcc*prcc))) - qt(0.975,df) ) 
}
psig <- uniroot(tval, c(0,1), N-5)$root ### df deve ser alterado NA MAO
abline(h=0);abline(h=psig,lty=2);abline(h=-psig,lty=2)

# TODO: che cazzo???
chisq.test(table(res,r))
chisq.test(res,Time)
chisq.test(res,k)

Ntest <- 5
cats <- qunif((1:Ntest)/Ntest, 0.25, 2)
rcats <- array(dim=N)
for (i in Ntest:1)
	rcats[r < cats[i]] <- i
rescats <-array(2,dim=N)
rescats[res<median(res)] <- 1
table(rescats,rcats)
chisq.test(table(rescats,rcats))
kruskal.test(res,rcats)

Kruskal <- array(dim=15)
for (Ntest in (2:15)) {
	cats <- qunif((1:Ntest)/Ntest, 1, 10)
	xcats <- array(dim=N)
	for (i in Ntest:1)
		xcats[x < cats[i]] <- i
	kruskal.test(res,xcats)$p.value -> Kruskal[Ntest]
}
Kruskal
plot(Kruskal, ylim=c(0,1))
abline(h=0.05)


Kruskal <- array(dim=15)
for (Ntest in (2:15)) {
	cats <- qunif((1:Ntest)/Ntest, 10, 50)
	kcats <- array(dim=N)
	for (i in Ntest:1)
		kcats[k < cats[i]] <- i
	kruskal.test(res,kcats)$p.value -> Kruskal[Ntest]
}
Kruskal
plot(Kruskal, ylim=c(0,1))
abline(h=0.05)

cor(res,k)
anova(lm(res~xcats))
