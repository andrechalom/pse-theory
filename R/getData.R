source("pse.R")
N <- 650
# Entradas da matriz:
# Vindos de Silva Matos, Freckleton & Watkinson 1999
dados <- read.csv("leslie.csv", header=TRUE)
# Media para cada entrada da matriz
means<-tapply(dados$VALUE, INDEX=list(dados$TO, dados$FROM), FUN=mean)
sds<-tapply(dados$VALUE, INDEX=list(dados$TO, dados$FROM), FUN=sd)

p11<- LHSsample (N, qprob=qnorm, means[1,1], sds[1,1])
# p17 tem uma variacao muito ampla, reduzindo arbitrariamente
p17<- LHSsample (N, qprob=qnorm, means[1,7], sds[1,7] * 0.7)
p21<- LHSsample (N, qprob=qnorm, means[2,1], sds[2,1])
p22<- LHSsample (N, qprob=qnorm, means[2,2], sds[2,2])
p32<- LHSsample (N, qprob=qnorm, means[3,2], sds[3,2])
p33<- LHSsample (N, qprob=qnorm, means[3,3], sds[3,3])
p43<- LHSsample (N, qprob=qnorm, means[4,3], sds[4,3])
p44<- LHSsample (N, qprob=qnorm, means[4,4], sds[4,4])
p54<- LHSsample (N, qprob=qnorm, means[5,4], sds[5,4])
p55<- LHSsample (N, qprob=qnorm, means[5,5], sds[5,5])
p65<- LHSsample (N, qprob=qnorm, means[6,5], sds[6,5])
p66<- LHSsample (N, qprob=qnorm, means[6,6], sds[6,6])
p76<- LHSsample (N, qprob=qnorm, means[7,6], sds[7,6])
# A media reportada para p77 eh considerada "alta demais" pelos autores.
# Correcao arbitraria
p77<- LHSsample (N, qprob=qnorm, means[7,7]-0.015, sds[7,7])

# Media de Gm retirada do paper, sd chutado (impossivel reconstruir)
Gm <- LHSsample (N, qprob=qnorm, 0.486, 0.15)
a  <- LHSsample (N, qprob=qnorm, 0.01228, 0.01228/3)
z  <- LHSsample (N, qprob=qnorm, 7, 1)
 
vars <- data.frame(p11,p17,p21,p22,p32,p33,p43,p44,p54,p55,p65,p66,p76,p77, Gm, a, z)
# Cortando dados impossiveis:
vars[vars<0.001] <- 0.001
vars[vars>0.999 & vars < 2] <- 0.999
vars <- LHScorcorr(vars)
save(N, vars, file="leslie.Rdata")
