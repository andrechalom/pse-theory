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
# Total de amostras sera 100:
N <- 100
# r eh uniformemente distribuido entre 0 e 2
# os outros parametros seguem logica semelhante
# (Samples para outras distribuicoes serao implementados posteriormente!)
r <- sampleunif(N, 0, 2, "r")
k <- sampleunif(N, 10, 50, "k")
x <- sampleunif(N, 1, 10, "x")
# O tempo final eh distribuido com forma normal de media = 100 e sd = 20
# TODO: samplenorm estah tomando os quantis, e nao uma amostra aleatoria!
# Como tomar uma amostra aleatoria de um intervalo de maneira eficiente??
Time <- samplenorm(N, 100, 20, "T")

# Correlacao entre as variaveis geradas. Veja que ela nao eh nec. 0
cor(cbind(r,k,Time,x))

# TODO: correcao de correlacoes!
# LHS(r,k,Time,x, <matriz de cor esperada>)

# Rodando o modelo para os diferentes parametros
# O objeto res contem a saida do modelo, no nosso caso, a populacao final
res <- mapply(modelRun, x, r, k, Time)

### Plots e analises
# TODO: Gerar uma estrutura mais generica para plotar graficos

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
# WORK IN PROGRESS!!!
# COMO EXTRAIR UMA "FATIA" MEANINGFUL DO CUBO??
predPlot(list(r,x,Time,k), res, modellm)

# Mais um modelo linear, agora restrito a regiao r > 1
modellm <- lm (res ~ r + x + Time + k, subset=(r>0.0001))
print (modellm)
print (summary(modellm))
predPlot(list(r,x,Time,k), res, modellm)

