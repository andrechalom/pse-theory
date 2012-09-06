source("pse.R")
N <- 650
# Entradas da matriz:
# Vindos de Silva Matos, Freckleton & Watkinson 1999
dados <- read.csv("leslie.csv", header=TRUE)
# Media para cada entrada da matriz
means<-tapply(dados$VALUE, INDEX=list(dados$TO, dados$FROM), FUN=mean)
sds<-tapply(dados$VALUE, INDEX=list(dados$TO, dados$FROM), FUN=sd)

P1<- LHSsample (N, qprob=qnorm, means[1,1], sds[1,1])
# p17 tem uma variacao muito ampla, reduzindo arbitrariamente
F7<- LHSsample (N, qprob=qnorm, means[1,7], sds[1,7] * 0.7)
G1<- LHSsample (N, qprob=qnorm, means[2,1], sds[2,1])
P2<- LHSsample (N, qprob=qnorm, means[2,2], sds[2,2])
G2<- LHSsample (N, qprob=qnorm, means[3,2], sds[3,2])
P3<- LHSsample (N, qprob=qnorm, means[3,3], sds[3,3])
G3<- LHSsample (N, qprob=qnorm, means[4,3], sds[4,3])
P4<- LHSsample (N, qprob=qnorm, means[4,4], sds[4,4])
G4<- LHSsample (N, qprob=qnorm, means[5,4], sds[5,4])
P5<- LHSsample (N, qprob=qnorm, means[5,5], sds[5,5])
G5<- LHSsample (N, qprob=qnorm, means[6,5], sds[6,5])
P6<- LHSsample (N, qprob=qnorm, means[6,6], sds[6,6])
G6<- LHSsample (N, qprob=qnorm, means[7,6], sds[7,6])
# A media reportada para p77 eh considerada "alta demais" pelos autores.
# Correcao arbitraria
P7<- LHSsample (N, qprob=qnorm, means[7,7]-0.015, sds[7,7])

# Media de Gm retirada do paper, sd chutado (impossivel reconstruir)
Gm <- LHSsample (N, qprob=qnorm, 0.486, 0.15)
a  <- LHSsample (N, qprob=qnorm, 0.01228, 0.01228/3)
z  <- LHSsample (N, qprob=qnorm, 7, 1)
 
vars <- data.frame(P1,F7,G1,P2,G2,P3,G3,P4,G4,P5,G5,P6,G6,P7, Gm, a, z)
# Cortando dados impossiveis:
vars[vars<0.001] <- 0.001
vars[vars>0.999 & vars < 2] <- 0.999
vars <- LHScorcorr(vars)
save(N, vars, file="leslie.Rdata")
