### Exemplos de utilizacao do pse:
# Primeiro, rodamos o modelo:
source("modelrun.R")

# Plot de sobrevivencia
# Nesta funcao, testamos uma condicao da saida do modelo: populacao
# final maior do que 0 (mais um errinho numerico)
# Note que as funcoes de plot aceitam parametros graficos com algumas
# condicoes (mais sobre isso no manual)
# A funcao plota + onde o teste eh verdadeiro, - onde eh falso
testPlot(list(r,x,Time,k), test=res>0.0001, mar=c(4,4,1,2))

testPlot(list(r,x,Time,k), test=res>0.0001, convex=TRUE, mar=c(4,4,1,2))

?deldir
?plot.triang.list
?plot.tile.list
