matriz = matrix(c(0, 2, 0.6, 0), ncol=2, byrow=TRUE)
lambda = eigen(matriz)$values[1]

# Verossimilhanca para gerar as	distribuicoes:
# primeiro censo, s=10, total=30
# segundo censo, s=7, total=25

# binomial com uma observacao, s = 10, total = 30
L = function (p) { return (931395465 * p^10 * (1-p)^20) }
curve(L(x), add=T)
curve (dbinom(10, 30, p=x))
D = function (p) dbinom(10, 30, p)
mult = integrate(D,0,1)$value
cD = function (p) integrate(D, 0, p)$value/mult
cD = Vectorize(cD)
curve(cD)

# veros. cumulativa da binomial:
cL = function (p) { return (integrate(L, 0, p)$value)}
cL = Vectorize(cL)
curve(cL(x))

# Usando dois censos independentes, com as hipoteses FORTES de que
# 1) a probabilidade de transicao se mantem constante a cada ano
# 2) a probabilidade de um evento com um dado individuo eh independente a cada ano

# Segundo censo
L2 = function (p) {return (12498200 * p^7 * (1-p)^18) } 
cL2 = function (p) { return (integrate(L2, 0, p)$value)}
cL2 = Vectorize(cL2)
curve(cL2(x), col="blue", add=T)

# Censos combinados
comboL = function(p) {L(p)*L2(p)/3.045803}
combo = function (p) {return (integrate(comboL,0,p)$value)}
combo = Vectorize(combo)
curve(combo, col="red", add=T)

# Ou podemos aproveitar que se X1 e X2 ~ binom, X1+X2 eh binom tb
somaL = function (p) {return (3821903815930200 * p^17 * (1-p)^(55-17)) }
soma = function (p) {return (integrate(somaL,0,p)$value)}
soma = Vectorize(soma)
curve(soma, col="orange", add=T)
# Vejam, eh igual =D
