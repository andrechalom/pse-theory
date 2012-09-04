source("pse.R")
library(sensitivity)
# Hipercubo gerado de tamanho 650
# Analises realizadas para modelo independente de densidade
load("Independent.Rdata")
# Analises realizadas para modelo dependente de densidade
load("Dependent.Rdata")

# Plots
par(mfrow=c(2,1), cex=0.8)
plot(imorr, xlim=c(0, 1.2*max(apply(imorr$ee, 2, function(x) mean(abs(x)))))); title("Morris screening, independent")
plot(dmorr, xlim=c(0, 1.2*max(apply(dmorr$ee, 2, function(x) mean(abs(x)))))); title("Morris screening, dependent")

par(mfrow=c(1,2))
plot(ecdf(ires), do.points=FALSE, main="ECDF, independent")
plot(ecdf(dres), do.points=FALSE, main="ECDF, dependent")

corPlot(ivars, ires)
corPlot(dvars, dres)

par(mfrow=c(2,1), cex=1, cex.axis=0.8)
plot(iprcc); abline(h=0, lty=3)
plot(dprcc); abline(h=0, lty=3)
# p77 e p11 pulam!

par(mfrow=c(2,1), cex.axis=0.8)
plot(ifast); title(main="eFAST")
plot(dfast); title(main="eFAST")
