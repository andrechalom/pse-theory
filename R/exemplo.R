### Exemplos de utilizacao do pse:
source("pse.R")
### Modelo numerico: crescimento logistico simples
# O modelo tem 4 parametros: r, K, populacao inicial e tempo final
# Vamos estudar qual a influencia de cada um deles
modelRun <- function (Xo, r, K, Time) {
	X <- Xo
	for (i in 0:Time) {
		X <- r*X[length(X)]*(1-X[length(X)]/K)
	}
	return (X)
}

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

testPlot(list(r,x,Time,k), test=res>0.0001, convex=TRUE, mar=c(4,4,1,2))
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

# O teste de KW recebe como input um vetor x e um vetor
# de fatores g. Para simplificar o uso do teste no caso
# de variaveis continuas, criamos uma funcao acessoria:
kruskal.test.c <- function (res, x, N) {
	cats <- cut(x, breaks=quantile(x, seq(0,1,1/N)))
	return(kruskal.test(res~cats))
}

# Teste de Kruskal-Wallis p/ variavel r:
kruskal.test.c(res,r,5)

# Teste de KW para N=2 ateh 15:
Kruskal <- array(dim=15)
for (Ntest in (2:15)) {
	kruskal.test.c(res,x,Ntest)$p.value -> Kruskal[Ntest]
}
Kruskal
plot(Kruskal, ylim=c(0,1))
abline(h=0.05) # Significancia usual


Kruskal <- array(dim=15)
for (Ntest in (2:15)) {
	kruskal.test.c(res,k,Ntest)$p.value -> Kruskal[Ntest]
}
Kruskal
plot(Kruskal, ylim=c(0,1))
abline(h=0.05)

cor(res,k)
anova(lm(res~x))


## ConvexHull
test
require(spatstat)
oneTestPlot(r,k, test=res>0.0001)
xTrue <- r[res<0.0001]; kk <- k[res<0.0001]
W <- owin(c(min(rr),max(rr)), c(min(kk),max(kk)))
P <- ppp(rr,kk,window=W)
C <-convexhull(P)
plot(C, add=T, col=cm.colors(1), density=10)
