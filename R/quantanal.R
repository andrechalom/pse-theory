### Exemplos de utilizacao do pse:
# Primeiro, rodamos o modelo:
source("modelrun.R")
# Analise de incerteza:
# Distribuicao de probabilidades de y = f(x)
hist(res[res>0.001])
plot(ecdf(res[res>0.0001]), do.points=FALSE)
par(mfrow=c(2,2))
# Aqui, N vale 20
plot(density(res[res>0.0001]))
# Entao, refinamos o numero de pontos:
source("asr.R")
# E vemos os graficos de res novamente:
plot(density(res[res>0.0001]))
source("asr.R")
plot(density(res[res>0.0001]))
source("asr.R")
# Aqui, N vale 540
plot(density(res[res>0.0001]))

# Alternativa: Load num sorteio jah realizado
# load("../Rdata/LHS1500.Rdata")

### Plots e analises
# Analise da correlacao entre o resultado e cada variavel:
# metodos de Pearson e Spearman e porcentagem da variancia de res
# explicada por efeitos lineares e monotonicos dos inputs
cor(cbind(res,r,k,Time,x))
cor(res,r)^2 + cor(res,k)^2 + cor(res,Time)^2 + cor(res,x)^2
cor(cbind(res,r,k,Time,x), method="spearman")
cor(res,r, method="spearman")^2 + cor(res,k, method="spearman")^2 + cor(res,Time, method="spearman")^2 + cor(res,x, method="spearman")^2

# TODO: Gerar uma estrutura mais generica para plotar graficos
# Plot da contribuicao de cada variavel, independentemente
# Os indicadores de significancia nao podem ser levados a serio demais...
corPlot(list(r,k,Time,x),res)

# Plot da populacao final
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

# PRCC: partial rank correlations
xmat <- cbind(r, k, Time, x)
colnames(xmat)<-c("r","k","Time","x")
PRCC <- pcor(res, xmat)
plot.pcor(PRCC, N, 4)

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

# Idem para k
Kruskal <- array(dim=15)
for (Ntest in (2:15)) {
	kruskal.test.c(res,k,Ntest)$p.value -> Kruskal[Ntest]
}
Kruskal
plot(Kruskal, ylim=c(0,1))
abline(h=0.05)
