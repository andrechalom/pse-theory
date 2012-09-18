source("pse.R")
library(msm)
N <- 650
# Entradas da matriz:
# Vindos de Silva Matos, Freckleton & Watkinson 1999
dados <- read.csv("leslie.csv", header=TRUE)
stasis <- dados[dados$TO == dados$FROM & dados$FROM != 7,]
growth <- dados[dados$TO == (dados$FROM +1),]

# Taxa de fecundidade
fertilitymean <- mean(dados$VALUE[dados$TO==1 & dados$FROM==7]) 
fertilitysd <- sd(dados$VALUE[dados$TO==1 & dados$FROM==7])
F7<- LHSsample (N, qprob=qtnorm, fertilitymean,fertilitysd, 0, 500)

# Calculamos as taxas de sobrevivencia:
sobrev <- data.frame(stasis$FROM, stasis$TO, stasis$YEAR, stasis$VALUE+growth$VALUE)
colnames(sobrev) <- c("FROM", "TO", "YEAR", "VALUE")
sobrevmeans<-diag(tapply(sobrev$VALUE, INDEX=list(sobrev$TO, sobrev$FROM), FUN=mean))
sobrevsds<-diag(tapply(sobrev$VALUE, INDEX=list(sobrev$TO, sobrev$FROM), FUN=sd))
# A media reportada para p77 eh considerada "alta demais" pelos autores. Correcao arbitraria
sobrevmeans[7] <- mean(dados$VALUE[dados$TO==7 & dados$FROM==7])-0.015 
sobrevsds[7] <- sd(dados$VALUE[dados$TO==7 & dados$FROM==7])
s1<- LHSsample (N, qprob=qtnorm, sobrevmeans[1], sobrevsds[1], 0, 1)
s2<- LHSsample (N, qprob=qtnorm, sobrevmeans[2], sobrevsds[2], 0, 1)
s3<- LHSsample (N, qprob=qtnorm, sobrevmeans[3], sobrevsds[3], 0, 1)
s4<- LHSsample (N, qprob=qtnorm, sobrevmeans[4], sobrevsds[4], 0, 1)
s5<- LHSsample (N, qprob=qtnorm, sobrevmeans[5], sobrevsds[5], 0, 1)
s6<- LHSsample (N, qprob=qtnorm, sobrevmeans[6], sobrevsds[6], 0, 1)
s7<- LHSsample (N, qprob=qtnorm, sobrevmeans[7], sobrevsds[7], 0, 1)

# Taxa de crescimento real:
realgrowth <- data.frame(stasis$FROM, stasis$TO, stasis$YEAR, (sobrev$VALUE-stasis$VALUE)/sobrev$VALUE)
colnames(realgrowth) <- c("FROM", "TO", "YEAR", "VALUE")
realgrowthmeans<-diag(tapply(realgrowth$VALUE, INDEX=list(realgrowth$TO, realgrowth$FROM), FUN=mean))
realgrowthsds<-diag(tapply(realgrowth$VALUE, INDEX=list(realgrowth$TO, realgrowth$FROM), FUN=sd))
g1<- LHSsample (N, qprob=qtnorm, realgrowthmeans[1], realgrowthsds[1], 0, 1)
g2<- LHSsample (N, qprob=qtnorm, realgrowthmeans[2], realgrowthsds[2], 0, 1)
g3<- LHSsample (N, qprob=qtnorm, realgrowthmeans[3], realgrowthsds[3], 0, 1)
g4<- LHSsample (N, qprob=qtnorm, realgrowthmeans[4], realgrowthsds[4], 0, 1)
g5<- LHSsample (N, qprob=qtnorm, realgrowthmeans[5], realgrowthsds[5], 0, 1)
g6<- LHSsample (N, qprob=qtnorm, realgrowthmeans[6], realgrowthsds[6], 0, 1)

# Media de gm retirada do paper, sd chutado (impossivel reconstruir)
gmmean <- 0.486/(0.486+mean(dados$VALUE[dados$TO==1 & dados$FROM==1])) 
gm <- LHSsample (N, qprob=qtnorm, gmmean, 0.06, 0, 1)
a  <- LHSsample (N, qprob=qtnorm, 0.01228, 0.01228/3, 0, 1)
z  <- LHSsample (N, qprob=qtnorm, 7, 1, 0, 20)

# Para rodar somente se necessario:
#vars <- data.frame(s1,F7,g1,s2,g2,s3,g3,s4,g4,s5,g5,s6,g6,s7, gm, a, z)
#vars <- LHScorcorr(vars)
#save(N, vars, file="leslie.Rdata")
